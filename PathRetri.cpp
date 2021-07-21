/*
 * PathRetri.cpp
 *
 *  Created on: 7 May 2021
 *      Author: zhangmengxuan
 */
#include "head.h"

/***CH Path Retrieval***/
vector<int> Graph::PathRetriCH(int ID1, int ID2){
	vector<int> path;
	int meetNodeID;
	if(ID1==ID2) return path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return path;

	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);
	//predecessor
	vector<int> PreForward(nodenum, -1);
	vector<int> PreBackward(nodenum, -1);
	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

	int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

	while(!fHeapForward.empty() || !fHeapBackward.empty() )
	{
		if(bF && bB)
			break;
		if(bF && fHeapBackward.empty())
			break;
		if(bB && fHeapForward.empty())
			break;
		//Forward Search
		if(!fHeapForward.empty() && !bF)
		{
			fHeapForward.extract_min(topNodeIDForward, topDisForward);

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
					meetNodeID=topNodeIDForward;
				}
			}

			for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second.first;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
						PreForward[neighborNodeID]=topNodeIDForward;
					}
				}
			}
		}

		//Backward Search
		if(!fHeapBackward.empty() && !bB)
		{
			fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

			if(vDistanceBackward[topNodeIDBackward] > d)
				bB = true;

			vVisitedB[topNodeIDBackward] = true;

			if(vVisitedF[topNodeIDBackward]){
				int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
				if(distTmp<d){
					d=distTmp;
					meetNodeID=topNodeIDBackward;
				}
			}

			for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
				neighborNodeID = (*in).first;
				neighborLength = (*in).second.first;

				int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
				if(!vVisitedB[neighborNodeID]){
					if(vDistanceBackward[neighborNodeID]>db){
						vDistanceBackward[neighborNodeID] = db;
						fHeapBackward.update(neighborNodeID, db);
						PreBackward[neighborNodeID]=topNodeIDBackward;
					}
				}
			}
		}
	}

	//firstly, get the concise path
	path.push_back(meetNodeID);
	if(meetNodeID!=ID1){
		if(PreForward[meetNodeID]!=-1){
			int preFor=PreForward[meetNodeID];
			while(preFor!=ID1){
				path.insert(path.begin(),preFor);
				preFor=PreForward[preFor];
			}
			path.insert(path.begin(),preFor);
		}else
			path.insert(path.begin(),ID1);
	}
	if(meetNodeID!=ID2){
		if(PreBackward[meetNodeID]!=-1){
			int preBac=PreBackward[meetNodeID];
			while(preBac!=ID2){
				path.push_back(preBac);
				preBac=PreBackward[preBac];
			}
			path.push_back(preBac);
		}else
			path.push_back(ID2);
	}
	//Then, get the whole detailed path
	for(int i=0;i<path.size()-1;i++){
		int currentID = path[i];
		while(OutEdgesM[currentID].find(path[i+1])!=OutEdgesM[currentID].end()){//find the original adjacent node
			int MnodeID=OutEdgesM[currentID][path[i+1]];
			path.insert(path.begin()+i+1, MnodeID);//add the MnodeID after the current i-th node
		}
	}
	//cout<<"path length "<<path.size()<<endl;
	return path;
}

void Graph::CHconsPath(string orderfile){
	ifstream IF(orderfile);
	if(!IF){
		cout<<"Cannot open Map "<<orderfile<<endl;
	}
	NodeOrder.assign(nodenum, -1);
	vNodeOrder.assign(nodenum, -1);
	int num, nodeID, nodeorder;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>nodeID>>nodeorder;
		NodeOrder[nodeID]=nodeorder;
		if(nodeorder!=-1){
			vNodeOrder[nodeorder]=nodeID;
		}
	}

	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);

	map<int, vector<int>> mi;
	SCconNodesMT.assign(nodenum, mi);
	map<int,int> mapinfor;//(endpoint, supportnode)
	OutEdgesM.assign(nodenum, mapinfor);

	//initialize E
	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
	}

	vector<bool> exist; exist.assign(nodenum,true);
	//vector<bool> change; change.assign(nodenum,false);

	//cout<<"Begin to contract"<<endl;
	for(int nodeorder=0;nodeorder<nodenum;nodeorder++){
		int x=vNodeOrder[nodeorder];
		if(x!=-1){//to identify and exclude the isolated vertices
			exist[x]=false;

			vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

			for(auto it=E[x].begin();it!=E[x].end();it++){
				if(exist[(*it).first]){
					Neigh.push_back(*it);
				}
			}
			NeighborCon[x].assign(Neigh.begin(),Neigh.end());

			for(int i=0;i<Neigh.size();i++){
				int y=Neigh[i].first;
				deleteEorder(x,y);
				//change[y]=true;
			}

			if(Neigh.size()<=100){
				//single thread
				for(int i=0;i<Neigh.size();i++){
					for(int j=i+1;j<Neigh.size();j++){
						bool support=insertEPath(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
						if(support){
							OutEdgesM[Neigh[i].first][Neigh[j].first]=x;
							OutEdgesM[Neigh[j].first][Neigh[i].first]=x;
						}
						/*if(Neigh[i].first<Neigh[j].first)
							SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
						else if(Neigh[j].first<Neigh[i].first)
							SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);*/
					}
				}
			}else{
				//multiple thread
				int step=Neigh.size()/threadnum;
				boost::thread_group thread;
				for(int i=0;i<threadnum;i++){
					pair<int,int> p;
					p.first=i*step;
					if(i==threadnum-1)
						p.second=Neigh.size();
					else
						p.second=(i+1)*step;
					thread.add_thread(new boost::thread(&Graph::NeighborComPath, this, boost::ref(Neigh), p, x));
				}
				thread.join_all();
			}

		}
	}
}

void Graph::NeighborComPath(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		w1=Neighvec[k].second.first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			w2=Neighvec[h].second.first;
			bool support=insertEMTPath(ID1, ID2, w1+w2);
			if(support){
				OutEdgesM[ID1][ID2]=x;
			}
			/*if(ID1<ID2)
				SCconNodesMT[ID1][ID2].push_back(x);*/
		}
	}
	sm->notify();
}

bool Graph::insertEPath(int u,int v,int w){
	bool support=false;
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
		support=true;
	}
	else{
		if(E[u][v].first>w){
			E[u][v]=make_pair(w,1);
			support=true;
		}
		else if(E[u][v].first==w){
			E[u][v].second+=1;
		}
	}

	if(E[v].find(u)==E[v].end()){
		E[v].insert(make_pair(u,make_pair(w,1)));
	}
	else{
		if(E[v][u].first>w)
			E[v][u]=make_pair(w,1);
		else if(E[v][u].first==w)
			E[v][u].second+=1;
	}
	return support;
}

bool Graph::insertEMTPath(int u,int v,int w){
	bool support=false;
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
		support=true;
	}
	else{
		if(E[u][v].first>w){
			E[u][v]=make_pair(w,1);
			support=true;
		}
		else if(E[u][v].first==w){
			E[u][v].second+=1;
		}

	}
	return support;
}

/**************************************************************/
/***H2H Path Retrieval***/
vector<int> Graph::PathRetriH2H(int ID1, int ID2){
	vector<int> Path;
	if(ID1==ID2) return Path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return Path;

	int r1=rank[ID1], r2=rank[ID2];
	int LCA=LCAQuery(r1,r2); //cout<<"LCA vertex "<<Tree[LCA].uniqueVertex<<endl;

	int tmp=INF;//shortest distance
	int LCAonSP;//the highest vertex on shortest path
	int seq=INF;//the ancestor position of the highest node on shortest path

	if(LCA==r1){
		//cout<<"case A"<<endl;
		Path=PathRetriH2Hpartial(ID2,ID1);
	}
	else if(LCA==r2){
		//cout<<"case B"<<endl;
		Path=PathRetriH2Hpartial(ID1,ID2);
	}
	else{
		//cout<<"case C"<<endl;
		for(int i=0;i<Tree[LCA].pos.size();i++){

			if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
				tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
				if(i==Tree[LCA].pos.size()-1)
					LCAonSP=Tree[LCA].uniqueVertex;
				else
					LCAonSP=Tree[LCA].vert[i].first;
				seq=Tree[LCA].pos[i];
			}else if(tmp==Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
				if(Tree[LCA].pos[i]<seq){
					if(i==Tree[LCA].pos.size()-1)
						LCAonSP=Tree[LCA].uniqueVertex;
					else
						LCAonSP=Tree[LCA].vert[i].first;
					seq=Tree[LCA].pos[i];
				}
			}

		}

		Path=PathRetriH2Hpartial(ID1,LCAonSP);

		vector<int> path1=PathRetriH2Hpartial(ID2, LCAonSP);
		for(int k=path1.size()-2;k>=0;k--){
			Path.push_back(path1[k]);
		}
	}
	cout<<"path length "<<Path.size()<<endl;
	return Path;
}

vector<int> Graph::PathRetriH2Hpartial(int ID1, int ID2){//ID2 is the ancestor of ID1
	vector<int> Path;
	if(ID1==ID2) return Path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return Path;

	int r1=rank[ID1], r2=rank[ID2];
	//int LCA=LCAQuery(r1,r2); //cout<<"Unique vertex of LCA "<<Tree[LCA].uniqueVertex<<endl;

	int tmp=INF;//shortest distance
	int LCAonSP;//the highest vertex on shortest path
	int seq=INF;//the ancestor position of the highest node on shortest path


		tmp = Tree[r1].dis[Tree[r2].pos.back()];
		seq=Tree[r2].pos.back();
		LCAonSP=ID2;
		int currentLCA=LCAonSP;

		int curseq;
		for(int k=0;k<Tree[rank[currentLCA]].pos.size()-1;k++){
			curseq=Tree[rank[currentLCA]].pos[k];
			if(Tree[r1].dis[curseq]+Tree[r2].dis[curseq]==tmp){
				if(curseq<seq){
					seq=curseq;
					LCAonSP=Tree[rank[currentLCA]].vert[k].first;
					//cout<<"////////////// "<<"LCA "<<LCAonSP<<endl;
				}
			}
		}

		//continue to find higher nodes
		while(LCAonSP!=currentLCA){
			currentLCA=LCAonSP;

			for(int k=0;k<Tree[rank[currentLCA]].pos.size()-1;k++){
				curseq=Tree[rank[currentLCA]].pos[k];
				if(Tree[r1].dis[curseq]+Tree[r2].dis[curseq]==tmp){
					if(curseq<seq){
						seq=curseq;
						LCAonSP=Tree[rank[currentLCA]].vert[k].first;
					}
				}
			}
			//cout<<"************* "<<"LCA "<<LCAonSP<<endl;
		}


	//get the concise path
	if(LCAonSP==ID2){
		int curID=ID1;
		while(Tree[rank[curID]].piv[seq]!=-1){
			Path.push_back(curID);
			curID=Tree[rank[curID]].piv[seq];
		}
		Path.push_back(curID);

		Path.push_back(LCAonSP);
	}else{
		int curID=ID1;
		while(Tree[rank[curID]].piv[seq]!=-1){
			Path.push_back(curID);
			curID=Tree[rank[curID]].piv[seq];
		}
		Path.push_back(curID);
		Path.push_back(LCAonSP);

		vector<int> path2;
		path2.clear();
		//from LCAonSP to ID2
		curID=ID2;
		while(Tree[rank[curID]].piv[seq]!=-1){
			path2.insert(path2.begin(), curID);
			curID=Tree[rank[curID]].piv[seq];
		}
		path2.insert(path2.begin(), curID);

		Path.insert(Path.end(), path2.begin(), path2.end());
	}

	//get the detailed path
	for(int i=0;i<Path.size()-1;i++){
		int currentID = Path[i];
		while(OutEdgesM[currentID].find(Path[i+1])!=OutEdgesM[currentID].end()){//find the original adjacent node
			int MnodeID=OutEdgesM[currentID][Path[i+1]];
			Path.insert(Path.begin()+i+1, MnodeID);//add the MnodeID after the current i-th node
		}
	}
	//cout<<"path length "<<Path.size()<<endl;
	return Path;
}

void Graph::H2HconPath(string orderfile){
	CHconsPath(orderfile);
	makeTree();
	makeIndexPath();
}

void Graph::makeIndexPath(){
	makeRMQ();

	//initialize
	vector<int> list; //list.clear();
	list.push_back(Tree[0].uniqueVertex);
	Tree[0].pos.clear();
	Tree[0].pos.push_back(0);

	for(int i=0;i<Tree[0].ch.size();i++){
		makeIndexDFSPath(Tree[0].ch[i],list);
	}
}

void Graph::makeIndexDFSPath(int p, vector<int>& list){
	//initialize
	int NeiNum=Tree[p].vert.size();
	Tree[p].pos.assign(NeiNum+1,0);
	Tree[p].dis.assign(list.size(),INF);
	Tree[p].piv.assign(list.size(),-1);//-1 means the super edge itself is the shortest distance
	Tree[p].cnt.assign(list.size(),0);
	Tree[p].FN.assign(list.size(),true);

	//pos
	//map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
	for(int i=0;i<NeiNum;i++){
		for(int j=0;j<list.size();j++){
			if(Tree[p].vert[i].first==list[j]){
				Tree[p].pos[i]=j;
				Tree[p].dis[j]=Tree[p].vert[i].second.first;
				Tree[p].cnt[j]=1;
				break;
			}
		}
	}
	Tree[p].pos[NeiNum]=list.size();


	//dis
	for(int i=0;i<NeiNum;i++){
		int x=Tree[p].vert[i].first;
		int disvb=Tree[p].vert[i].second.first;
		int k=Tree[p].pos[i];//the kth ancestor is x

		for(int j=0;j<list.size();j++){
			int y=list[j];//the jth ancestor is y

			int z;//the distance from x to y
			if(k!=j){
				if(k<j)
					z=Tree[rank[y]].dis[k];
				else if(k>j)
					z=Tree[rank[x]].dis[j];

				if(Tree[p].dis[j]>z+disvb){
					Tree[p].dis[j]=z+disvb;
					Tree[p].FN[j]=false;
					Tree[p].cnt[j]=1;
					Tree[p].piv[j]=x;
				}else if(Tree[p].dis[j]==z+disvb){
					Tree[p].cnt[j]+=1;
				}
			}
		}
	}

	//nested loop
	list.push_back(Tree[p].uniqueVertex);
	for(int i=0;i<Tree[p].ch.size();i++){
		makeIndexDFSPath(Tree[p].ch[i],list);
	}
	list.pop_back();
}

/**********************************************************************/

/****PLL Path Retrieval****/

void Graph::PLLconsPath(string orderfile, string indexfile){
	ifstream IF(orderfile);
	if(!IF){
		cout<<"Cannot open Map "<<orderfile<<endl;
	}
	NodeOrder.assign(nodenum, 0);
	vNodeOrder.assign(nodenum, 0);
	int num, nodeID, nodeorder;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>nodeID>>nodeorder;
		NodeOrder[nodeID]=nodeorder;
		vNodeOrder[nodeorder]=nodeID;
	}

	unordered_map<int,int> m;
	labelPos.assign(nodenum,m);
	vector<pair<int,int>> m1;
	label.assign(nodenum,m1);
	vector<pair<int,pair<int,int>>> m2;
	labelPre.assign(nodenum, m2);

	fstream file;
	file.open(indexfile);
	if(!file){
		int ID, nodeid, dist, prece;
		int cnt=0;
		for(int i=nodenum-1;i>=0;i--){
			ID=vNodeOrder[i];
			vector<pair<int,pair<int,int>>> vp;
			DijkPLLPath(ID,vp);
			for(int j=0;j<vp.size();j++){
				nodeid=vp[j].first;
				dist=vp[j].second.first;
				prece=vp[j].second.second;
				label[nodeid].push_back(make_pair(ID, dist));
				labelPre[nodeid].push_back(make_pair(ID,make_pair(dist,prece)));
				labelPos[nodeid][ID]=label[nodeid].size()-1;
			}
			cnt+=1;
			if(cnt%10000==0) cout<<"cnt "<<cnt<<endl;
		}
		PLLPathIndexWrite(indexfile);
		PLLPathIndexRead(indexfile);
	}else
		PLLPathIndexRead(indexfile);
}

void Graph::DijkPLLPath(int nodeID, vector<pair<int,pair<int,int>>> &vp){
	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);
	vector<int> prece(nodenum, -1);
	int topNodeID, topWeight;

	vector<int> singlelabel(nodenum, INF);
	for(int i=0;i<label[nodeID].size();i++)
		singlelabel[label[nodeID][i].first]=label[nodeID][i].second;

	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(nodeID,0);
	distance[nodeID]=0;
	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topWeight);
		closed[topNodeID]=true;

		int TempDis=INF;
		for(int j=0;j<label[topNodeID].size();j++){
			int hubID=label[topNodeID][j].first;
			int dis=label[topNodeID][j].second;
			if(singlelabel[hubID]>topWeight || dis>topWeight)
				continue;
			else{
				int disnew=dis+singlelabel[hubID];
				if(disnew<TempDis)
					TempDis=disnew;
			}
		}

		if(topWeight>=TempDis)
			continue;

		if(prece[topNodeID]==nodeID)
			vp.push_back(make_pair(topNodeID, make_pair(topWeight,-1)));
		else
			vp.push_back(make_pair(topNodeID, make_pair(topWeight, prece[topNodeID])));

		int NNodeID, NWeight;
		for(int k=0;k<Neighbor[topNodeID].size();k++){
			NNodeID=Neighbor[topNodeID][k].first;
			NWeight=topWeight+Neighbor[topNodeID][k].second;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeight){
					distance[NNodeID]=NWeight;
					pqueue.update(NNodeID, NWeight);
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}
}

void Graph::PLLPathIndexWrite(string file){
	cout<<"Index Writing"<<endl;
	ofstream outfile(file);

	vector<pair<int,pair<int,int>>>::iterator ivp;
	for(int i=0;i<nodenum;i++){
		outfile<<i<<"\t"<<label[i].size();
		for(ivp=labelPre[i].begin();ivp!=labelPre[i].end();ivp++)
			outfile<<"\t"<<(*ivp).first<<"\t"<<(*ivp).second.first<<"\t"<<(*ivp).second.second;
		outfile<<endl;
	}
	outfile.close();
}

void Graph::PLLPathIndexRead(string file){
	cout<<"Index Reading"<<endl;
	ifstream infile(file);

	int nodeid,hopnum,hop,d,pre;
	for(int i=0;i<nodenum;i++){
		infile>>nodeid>>hopnum;
		for(int j=0;j<hopnum;j++){
			infile>>hop>>d>>pre;
			label[nodeid].push_back(make_pair(hop,d));
			labelPre[nodeid].push_back(make_pair(hop,make_pair(d,pre)));
			labelPos[nodeid][hop]=j;
		}
	}
}

vector<int> Graph::PathRetriPLL(int ID1, int ID2){
	vector<int> path;
	if(ID1==ID2) return path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return path;
	path.push_back(ID1);
	path.push_back(ID2);
	int distance=DistanceCompute(ID1,ID2);

	if(distance!=INF){
		for(int j=0;j<path.size()-1;j++){
			if(AdjacentNodes[path[j]].find(path[j+1])!=AdjacentNodes[path[j]].end())
				continue;
			while(AdjacentNodes[path[j]].find(path[j+1])==AdjacentNodes[path[j]].end()){
				vector<int> Points=PointCompute(path[j],path[j+1]);
				path.insert(path.begin()+j+1,Points.begin(),Points.end());
			}
		}
	}
	//cout<<"path retrieved with size "<<path.size()<<endl;
	return path;
}

int Graph::DistanceCompute(int ID1, int ID2){
	int shortestdis=INF;
	int size1=label[ID1].size();
	int size2=label[ID2].size();

	if(size1==0 && size2==0) return shortestdis;

	if(labelPos[ID1].find(ID2)!=labelPos[ID1].end()){
		return labelPre[ID1][labelPos[ID1][ID2]].second.first;
	}else if(labelPos[ID2].find(ID1)!=labelPos[ID2].end()){
		return labelPre[ID2][labelPos[ID2][ID1]].second.first;
	}else{
		if(size1<size2){
			for(int i=0;i<label[ID1].size();i++){
				int id=label[ID1][i].first;
				if(labelPos[ID2].find(id)!=labelPos[ID2].end()){
					int dis=label[ID1][i].second+label[ID2][labelPos[ID2][id]].second;
					if(dis<shortestdis){
						shortestdis=dis;
					}
				}
			}
		}else if(size2<size1){
			for(int j=0;j<label[ID2].size();j++){
				int id=label[ID2][j].first;
				if(labelPos[ID1].find(id)!=labelPos[ID1].end()){
					int dis=label[ID2][j].second+label[ID1][labelPos[ID1][id]].second;
					if(dis<shortestdis){
						shortestdis=dis;
					}
				}
			}
		}
		return shortestdis;
	}
}

vector<int> Graph::PointCompute(int ID1, int ID2){
	vector<int> points;
	int size1=label[ID1].size();
	int size2=label[ID2].size();

	if(size1==0 && size2==0) return points;

	if(labelPos[ID1].find(ID2)!=labelPos[ID1].end()){
		if(labelPre[ID1][labelPos[ID1][ID2]].second.second!=-1)
			points.push_back(labelPre[ID1][labelPos[ID1][ID2]].second.second);
		return points;
	}else if(labelPos[ID2].find(ID1)!=labelPos[ID2].end()){
		if(labelPre[ID2][labelPos[ID2][ID1]].second.second!=-1)
			points.push_back(labelPre[ID2][labelPos[ID2][ID1]].second.second);
		return points;
	}else{
		int shortestdis=INF; int hopID;
		if(size1<size2){
			for(int i=0;i<label[ID1].size();i++){
				int id=label[ID1][i].first;
				if(labelPos[ID2].find(id)!=labelPos[ID2].end()){
					int dis=label[ID1][i].second+label[ID2][labelPos[ID2][id]].second;
					if(dis<shortestdis){
						shortestdis=dis;
						hopID=id;
					}
				}
			}
		}else if(size2<size1){
			for(int j=0;j<label[ID2].size();j++){
				int id=label[ID2][j].first;
				if(labelPos[ID1].find(id)!=labelPos[ID1].end()){
					int dis=label[ID2][j].second+label[ID1][labelPos[ID1][id]].second;
					if(dis<shortestdis){
						shortestdis=dis;
						hopID=id;
					}
				}
			}
		}
		int id1=labelPre[ID1][labelPos[ID1][hopID]].second.second;
		if(id1!=-1) points.push_back(id1);
		points.push_back(hopID);
		int id2=labelPre[ID2][labelPos[ID2][hopID]].second.second;
		if(id2!=-1) points.push_back(id2);
		return points;
	}
}

/****************************************************************/
/***CHP path retrieval***/
vector<int>	Graph::PathRetriCHP(int ID1, int ID2){
	vector<int> path;
	if(ID1==ID2) return path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return path;
	int d=INF;
	int meetNodeID;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

	//closed or not
	vector<bool> vVisitedF(nodenum, false);
	vector<bool> vVisitedB(nodenum, false);
	//the existing shortest distance
	vector<int>	vDistanceForward(nodenum, INF);
	vector<int>	vDistanceBackward(nodenum, INF);
	//predecessor
	vector<int> PreForward(nodenum, -1);
	vector<int> PreBackward(nodenum, -1);
	//stop search or not
	bool bF = false;
	bool bB = false;
	vDistanceForward[ID1] = 0;
	vDistanceBackward[ID2] = 0;
	fHeapForward.update(ID1,0);
	fHeapBackward.update(ID2,0);

	int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

	while(!fHeapForward.empty() || !fHeapBackward.empty() )
	{
		if(bF && bB)
			break;
		if(bF && fHeapBackward.empty())
			break;
		if(bB && fHeapForward.empty())
			break;
		//Forward Search
		if(!fHeapForward.empty() && !bF)
		{
			fHeapForward.extract_min(topNodeIDForward, topDisForward);

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
					meetNodeID=topNodeIDForward;
				}
			}

			for(auto out=AdjaShort[topNodeIDForward].begin();out!=AdjaShort[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
						PreForward[neighborNodeID]=topNodeIDForward;
					}
				}
			}
		}

		//Backward Search
		if(!fHeapBackward.empty() && !bB)
		{
			fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

			if(vDistanceBackward[topNodeIDBackward] > d)
				bB = true;

			vVisitedB[topNodeIDBackward] = true;

			if(vVisitedF[topNodeIDBackward]){
				int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
				if(distTmp<d){
					d=distTmp;
					meetNodeID=topNodeIDBackward;
				}
			}

			for(auto in=AdjaShort[topNodeIDBackward].begin();in!=AdjaShort[topNodeIDBackward].end();in++){
				neighborNodeID = (*in).first;
				neighborLength = (*in).second;

				int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
				if(!vVisitedB[neighborNodeID]){
					if(vDistanceBackward[neighborNodeID]>db){
						vDistanceBackward[neighborNodeID] = db;
						fHeapBackward.update(neighborNodeID, db);
						PreBackward[neighborNodeID]=topNodeIDBackward;
					}
				}
			}
		}
	}

	//firstly, get the concise path
	path.push_back(meetNodeID);
	if(meetNodeID!=ID1){
		if(PreForward[meetNodeID]!=-1){
			int preFor=PreForward[meetNodeID];
			while(preFor!=ID1){
				path.insert(path.begin(),preFor);
				preFor=PreForward[preFor];
			}
			path.insert(path.begin(),preFor);
		}else
			path.insert(path.begin(),ID1);
	}
	if(meetNodeID!=ID2){
		if(PreBackward[meetNodeID]!=-1){
			int preBac=PreBackward[meetNodeID];
			while(preBac!=ID2){
				path.push_back(preBac);
				preBac=PreBackward[preBac];
			}
			path.push_back(preBac);
		}else
			path.push_back(ID2);
	}
	//Then, get the whole detailed path
	for(int i=0;i<path.size()-1;i++){
		int currentID = path[i];
		while(OutEdgesM[currentID].find(path[i+1])!=OutEdgesM[currentID].end()){//find the original adjacent node
			int MnodeID=OutEdgesM[currentID][path[i+1]];
			path.insert(path.begin()+i+1, MnodeID);//add the MnodeID after the current i-th node
		}
	}

	return path;
}

int Graph::StartCHPPath(string orderfile){
	ifstream IF(orderfile);
	if(!IF){
		cout<<"Cannot open Map "<<orderfile<<endl;
	}
	NodeOrder.assign(nodenum, -1);
	vNodeOrder.assign(nodenum, -1);
	int num, nodeID, nodeorder;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>nodeID>>nodeorder;
		NodeOrder[nodeID]=nodeorder;
		if(nodeorder!=-1){
			vNodeOrder[nodeorder]=nodeID;
		}
	}

	vector<pair<int,int>> vecp;
	AdjaShort.assign(nodenum, vecp);
	AdjaShortR.assign(nodenum, vecp);
	vvNode.assign(Neighbor.begin(), Neighbor.end());
	map<int, vector<int>> vecm;
	SupportNodes.assign(nodenum, vecm);
	map<int,int> mapinfor;
	OutEdgesM.assign(nodenum, mapinfor);
	vector<pair<int,int>> nodeinfor;
	vvpShortCut.assign(nodenum, nodeinfor);
	vector<bool> vbVisited(nodenum, false);
	map<int, int> mPosition;
	mPosition[0]=nodenum-1;
	bool bUpdated;
	vector<pair<int,int>> vU,vW;

	int v;// current contracting vertex
	int NOwCount=0;
	for(int i=0;i<vNodeOrder.size()-1;i++){
		v=vNodeOrder[i];
		if(v!=-1){
			vU=vvNode[v];
			for(auto ivpr=vU.begin();ivpr!=vU.end();ivpr++){
				CHcontractionorderPath((*ivpr).first, v, vbVisited, (*ivpr).second, vU);
			}

			//deleting v from G
			for(auto ivp = vvNode[v].begin(); ivp != vvNode[v].end(); ivp++)
			{
				for(auto ivpr = vvNode[(*ivp).first].begin(); ivpr != vvNode[(*ivp).first].end(); ivpr++)
					if((*ivpr).first == v)
					{
						vvNode[(*ivp).first].erase(ivpr);
						break;
					}
			}
			vbVisited[v]=true;
		}
	}

	return 0;
}

int Graph::CHcontractionorderPath(int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW){
	benchmark::heap<2,int,int> Heap(nodenum);
	vector<int> vDistance(nodenum, INF);
	int topNodeID, topDistance, neighborNodeID, neighborLength;
	int w;
	vDistance[ID1]=0;
	Heap.update(ID1,0);
	int count=0;
	map<int,int> mWDistance;
	vector<pair<int,int>> vpWDistance;
	map<int,int> mDistance;
	int maxWDistance=-1;

	for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
		w=(*ivp).first;
		if(NodeOrder[w]>NodeOrder[ID1]){//to get rid of redundant computation
			if(vbVisited[w])
				continue;
			int d=(*ivp).second+dUV;
			mWDistance[w]=d;
			if(d>maxWDistance)
				maxWDistance=d;
			mDistance[w]=INF;
			//maintain the support node information
			if(w<ID1)
				SupportNodes[w][ID1].push_back(ID2);
			else
				SupportNodes[ID1][w].push_back(ID2);
		}
	}
	if(mDistance.empty()){
		//NOwCount+=1;
		return 0;
	}
	for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
		vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

	int dThreshold=maxWDistance;

	while(!Heap.empty()){
		Heap.extract_min(topNodeID, topDistance);
		if(vbVisited[topNodeID])
			continue;
		if(topDistance>dThreshold)
			break;
		for(auto ivp=vvNode[topNodeID].begin();ivp!=vvNode[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(vbVisited[neighborNodeID] || neighborNodeID==ID2)
				continue;
			int d = vDistance[topNodeID] + neighborLength;
			if(vDistance[neighborNodeID] ==INF && neighborNodeID != ID1)
			{
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}else if(vDistance[neighborNodeID]>d){
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
			else if(vDistance[neighborNodeID] == d)
			{
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}
		}
	}

	for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
	{
		if((*imDistance).second > mWDistance[(*imDistance).first])//add shortcuts
		{
			int w = (*imDistance).first;
			int distance = mWDistance[(*imDistance).first];

			vvpShortCut[ID1].push_back(make_pair(w, distance));
			vvNode[ID1].push_back(make_pair(w, distance));
			vvNode[w].push_back(make_pair(ID1, distance));
			OutEdgesM[ID1].insert(make_pair(w,ID2));
			OutEdgesM[w].insert(make_pair(ID1,ID2));
		}
	}

	return 0;
}

/****************************************************************/
/*Bi-Dijkstra path retrieval*/
vector<int> Graph::BiDijPath(int ID1, int ID2){
	vector<int> path;
	if(ID1==ID2) return path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return path;//to avoid the incorrectness caused by the isolated vertex
	benchmark::heap<2, int, int> queueF(nodenum), queueB(nodenum);
	queueF.update(ID1,0);
	queueB.update(ID2,0);

	vector<bool> closedF(nodenum, false), closedB(nodenum, false);
	vector<int> distanceF(nodenum, INF), distanceB(nodenum, INF);
	vector<int> PreF(nodenum, -1), PreB(nodenum, -1);

	distanceF[ID1]=0;
	distanceB[ID2]=0;
	int topNodeIDF, topNodeDisF, topNodeIDB, topNodeDisB;
	int NNodeIDF,NWeighF, NNodeIDB, NWeighB;

	int d=INF;//initialize d to infinite for the unreachable case
	int meetID;

	while(!queueF.empty() || !queueB.empty()){
		if(queueF.top()+queueB.top()>=d){
			break;
		}
		//forward
		queueF.extract_min(topNodeIDF, topNodeDisF);
		closedF[topNodeIDF]=true;
		for(auto itF=Neighbor[topNodeIDF].begin();itF!=Neighbor[topNodeIDF].end();itF++){
			NNodeIDF=(*itF).first;
			NWeighF=(*itF).second+topNodeDisF;
			if(closedB[NNodeIDF] && NWeighF+distanceB[NNodeIDF]<d){
				d=NWeighF+distanceB[NNodeIDF];
				meetID=topNodeIDF;
			}
			if(!closedF[NNodeIDF]){
				if(distanceF[NNodeIDF]>NWeighF){
					distanceF[NNodeIDF]=NWeighF;
					queueF.update(NNodeIDF, NWeighF);
					PreF[NNodeIDF]=topNodeIDF;
				}
			}
		}
		//backward
		queueB.extract_min(topNodeIDB, topNodeDisB);
		closedB[topNodeIDB]=true;
		for(auto itB=Neighbor[topNodeIDB].begin();itB!=Neighbor[topNodeIDB].end();itB++){
			NNodeIDB=(*itB).first;
			NWeighB=(*itB).second+topNodeDisB;
			if(closedF[NNodeIDB] && NWeighB+distanceF[NNodeIDB]<d){
				d=NWeighB+distanceF[NNodeIDB];
				meetID=topNodeIDB;
			}
			if(!closedB[NNodeIDB]){
				if(distanceB[NNodeIDB]>NWeighB){
					distanceB[NNodeIDB]=NWeighB;
					queueB.update(NNodeIDB, NWeighB);
					PreB[NNodeIDB]=topNodeIDB;
				}
			}
		}
	}

	path.push_back(meetID);

	if(PreF[meetID]!=-1){
		int preFor=PreF[meetID];
		while(preFor!=ID1){
			path.insert(path.begin(),preFor);
			preFor=PreF[preFor];
		}
		path.insert(path.begin(),preFor);
	}else
		path.insert(path.begin(),ID1);

	if(PreB[meetID]!=-1){
		int preBac=PreB[meetID];
		while(preBac!=ID2){
			path.push_back(preBac);
			preBac=PreB[preBac];
		}
		path.push_back(preBac);
	}else
		path.push_back(ID2);

	//cout<<"Path length "<<path.size()<<endl;
	return path;
}

/****************************************************************/
/*A* path retrieval*/
vector<int> Graph::AstarPath(int ID1, int ID2){
	vector<int> path;
	if(ID1==ID2) return path;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return path;
	benchmark::heap<2, int, int> pqueue(nodenum);
	int heurisDis=EuclideanDis(ID1,ID2);
	pqueue.update(ID1,heurisDis);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);//actual distance from source
	vector<int> pre(nodenum, -1);

	distance[ID1]=0;
	int topNodeID, topNodeDis;
	int NNodeID,NWeigh;

	int d=INF;//initialize d to infinite for the unreachable case

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		if(topNodeID==ID2){
			d=distance[ID2];
			break;
		}
		closed[topNodeID]=true;

		for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+distance[topNodeID];
			heurisDis=EuclideanDis(NNodeID,ID2);
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh+heurisDis);
					pre[NNodeID]=topNodeID;
				}
			}
		}
	}

	path.push_back(ID2);
	int preFor=pre[ID2];
	while(preFor!=ID1){
		path.insert(path.begin(),preFor);
		preFor=pre[preFor];
	}
	path.insert(path.begin(),ID1);
	return path;
}
