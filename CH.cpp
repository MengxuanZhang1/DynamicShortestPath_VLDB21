/*
 * CHcon.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "head.h"

vector<int> _DD,_DD2;
struct DegComp{
	int x;
	DegComp(int _x){
		x=_x;
	}
	bool operator< (const DegComp d) const{
		if(_DD[x]!=_DD[d.x])
			return _DD[x]<_DD[d.x];
		if(_DD2[x]!=_DD2[x])
			return _DD2[x]<_DD2[d.x];
		return x<d.x;
	}
};

void Graph::CHconsMTOrderGenerate(){
	int Twidth=0;//tree width
	//initialize SCconNodesMT
	map<int, vector<int>> mi;
	SCconNodesMT.assign(nodenum, mi);

	//initialize E
	map<int,pair<int,int>> m;
	E.assign(nodenum,m);
	for(int i=0;i<Neighbor.size();i++){
		for(int j=0;j<Neighbor[i].size();j++)
			E[i].insert(make_pair(Neighbor[i][j].first,make_pair(0,1)));
	}

	_DD.assign(nodenum,0);_DD2.assign(nodenum,0);
	DD.assign(nodenum,0); DD2.assign(nodenum,0);

	set<DegComp> Deg;
	int degree;
	for(int i=0;i<nodenum;i++){
		degree=Neighbor[i].size();
		if(degree!=0){
			_DD[i]=degree;
			_DD2[i]=degree;
			DD[i]=degree;
			DD2[i]=degree;
			Deg.insert(DegComp(i));
		}
	}

	vector<bool> exist; exist.assign(nodenum,true);
	vector<bool> change; change.assign(nodenum,false);

	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect); //NeighborCon.clear();
	//SCconNodes.clear();

	//cout<<"Begin to contract"<<endl;
	int count=0;

	while(!Deg.empty()){
		if(count%10000==0)
			cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
		count+=1;
		int x=(*Deg.begin()).x;

		while(true){
			if(change[x]){
				Deg.erase(DegComp(x));
				_DD[x]=DD[x];
				_DD2[x]=DD2[x];
				Deg.insert(DegComp(x));
				change[x]=false;
				x=(*Deg.begin()).x;
			}else
				break;
		}

		vNodeOrder.push_back(x);
		Deg.erase(Deg.begin());
		exist[x]=false;

		vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

		for(auto it=E[x].begin();it!=E[x].end();it++){
			if(exist[(*it).first]){
				Neigh.push_back(*it);
			}
		}

		if(Neigh.size()>Twidth)
			Twidth=Neigh.size();

		NeighborCon[x].assign(Neigh.begin(),Neigh.end());

		//multi threads for n^2 combination
		for(int i=0;i<Neigh.size();i++){
			int y=Neigh[i].first;
			deleteE(x,y);
			change[y]=true;
		}

		int stepf=Neigh.size()/threadnum;
		boost::thread_group threadf;
		for(int i=0;i<threadnum;i++){
			pair<int,int> p;
			p.first=i*stepf;
			if(i==threadnum-1)
				p.second=Neigh.size();
			else
				p.second=(i+1)*stepf;
			threadf.add_thread(new boost::thread(&Graph::NeighborComOrderGenerate, this, boost::ref(Neigh), p, x));
		}
		threadf.join_all();
	}

	NodeOrder.assign(nodenum,-1);
	for(int k=0;k<vNodeOrder.size();k++){
		NodeOrder[vNodeOrder[k]]=k;
	}
	cout<<"Finish Contract"<<" , treewidth "<<Twidth<<endl;
}


int Graph::StartCHOrderMT(string indexfile, string orderfile){
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	fstream file;
	file.open(indexfile);
	if(!file){
		t1=std::chrono::high_resolution_clock::now();
		CHconsorderMT(orderfile);
		t2=std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
		runT= time_span.count();
		cout<<"CH order MT "<<runT<<endl;
		writeShortCutCH(indexfile);
		ReadShortCutCH(indexfile);
	}
	else
	{
		ReadShortCutCH(indexfile);
	}

	return 0;
}

void Graph::insertEMTorder(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
	}
	else{
		if(E[u][v].first>w)
			E[u][v]=make_pair(w,1);
		else if(E[u][v].first==w)
			E[u][v].second+=1;
	}
}

void Graph::CHconsorderMT(string orderfile){
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
						insertEorder(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
						if(Neigh[i].first<Neigh[j].first)
							SCconNodesMT[Neigh[i].first][Neigh[j].first].push_back(x);//no direction
						else if(Neigh[j].first<Neigh[i].first)
							SCconNodesMT[Neigh[j].first][Neigh[i].first].push_back(x);
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
					thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
				}
				thread.join_all();
			}

		}
	}
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		w1=Neighvec[k].second.first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			w2=Neighvec[h].second.first;
			insertEMTorder(ID1, ID2, w1+w2);
			if(ID1<ID2)
				SCconNodesMT[ID1][ID2].push_back(x);
		}
	}
	sm->notify();
}

//the multiple thread of CH with pruning, order decided beforehand
int Graph::StartCHPOrderMT(string indexfile, string orderfile){
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

	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	fstream file;
	file.open(indexfile);
	if(!file){
		t1=std::chrono::high_resolution_clock::now();
		CHPconsorderMT();
		t2=std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
		runT= time_span.count();
		cout<<"CHP MT "<<runT<<endl;
		writeShortCutorder(indexfile);
		ReadShortCut(indexfile);
	}
	else
	{
		ReadShortCut(indexfile);
	}

	return 0;
}

void Graph::CHPconsorderMT(){
	vector<pair<int,int>> nodeinfor;
	vvpShortCut.assign(nodenum, nodeinfor);
	vector<bool> vbVisited(nodenum, false);
	map<int, int> mPosition;
	mPosition[0]=nodenum-1;
	bool bUpdated;
	vector<pair<int,int>> vU,vW;

	//cout<<"Building static CH "<<endl;

	int v;// current contracting vertex
	for(int i=0;i<vNodeOrder.size()-1;i++){
		//if(i%10000==0) cout<<"i=/////////////////////// "<<i<<endl;
		v=vNodeOrder[i];
		if(v!=-1){
			vU=vvNode[v];
			//cout<<"contracting node "<<v<<endl;


			//multiple thread
			vector<vector<pair<int,int>>> vvpResult;
			vector<pair<int,int>> vec;
			vec.clear();
			vvpResult.assign(vU.size(), vec);
			//filter the neighbors no need for contraction
			vector<vector<pair<int,int>>> vUeach;
			vUeach.assign(vU.size(), vec);
			int w;
			for(int l=0;l<vU.size();l++){
				int ID1=vU[l].first;
				for(auto ivp=vU.begin();ivp!=vU.end();ivp++){
					w=(*ivp).first;
					if(NodeOrder[w]>NodeOrder[ID1]){//to get rid of redundant computation
						if(vbVisited[w])
							continue;
						vUeach[l].push_back(*ivp);
						if(w<ID1)
							SupportNodes[w][ID1].push_back(v);
						else
							SupportNodes[ID1][w].push_back(v);
					}
				}
			}

			//better not single thread
			/*if(vU.size()<10){
				//single thread
				for(auto ivpr=vU.begin();ivpr!=vU.end();ivpr++){
					CHcontractionorder((*ivpr).first, v, vbVisited, (*ivpr).second, vU);
				}
			}else{*/
				boost::thread_group threadf;
				for(int k=0;k<vU.size();k++){
					if(vUeach[k].size()>0)
						threadf.add_thread(new boost::thread(&Graph::CHcontractionorderMT, this, k, vU[k].first, v, boost::ref(vbVisited), vU[k].second, boost::ref(vUeach[k]), boost::ref(vvpResult)));
				}
				threadf.join_all();
				//merge to new shortcut and write it back
				for(int k=0;k<vU.size();k++)
				{
					int	ID1 = vU[k].first;
					for(int j=0;j<vvpResult[k].size();j++){
						int w=vvpResult[k][j].first;
						int distance=vvpResult[k][j].second;
						vvpShortCut[ID1].push_back(make_pair(w, distance));
						vvNode[ID1].push_back(make_pair(w, distance));
						vvNode[w].push_back(make_pair(ID1, distance));
					}
				}
			//}

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

}

int Graph::CHcontractionorderMT(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult){
	benchmark::heap<2,int,int> Heap(nodenum);
	vector<int> vDistance(nodenum, INF);
	int topNodeID, topDistance, neighborNodeID, neighborLength;
	int w;
	vDistance[ID1]=0;
	Heap.update(ID1,0);
	map<int,int> mWDistance;
	vector<pair<int,int>> vpWDistance;
	map<int,int> mDistance;
	int maxWDistance=-1;

	for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
		w=(*ivp).first;
		//if(NodeOrder[w]>NodeOrder[ID1]){//to get rid of redundant computation
			//if(vbVisited[w])
				//continue;
			int d=(*ivp).second+dUV;
			mWDistance[w]=d;
			if(d>maxWDistance)
				maxWDistance=d;
			mDistance[w]=INF;
		//}
	}

	if(mDistance.empty()){
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
		if((*imDistance).second > mWDistance[(*imDistance).first])
		{
			int w = (*imDistance).first;
			int distance = mWDistance[(*imDistance).first];

			vvpResult[k].push_back(make_pair(w, distance));
		}
	}

	return 0;
}

int Graph::StartCHPOrderMTWP(string indexfile, string orderfile){
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
	unordered_map<pair<int,int>, int, hash_pair> pathinfvec;
	PathInfor.assign(nodenum, pathinfvec);
	vector<tri> vec;
	EdgeOnPath.assign(Edge.size(),vec);

	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	fstream file;
	file.open(indexfile);
	if(!file){
		t1=std::chrono::high_resolution_clock::now();
		CHPconsorderMTWP();
		t2=std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
		runT= time_span.count();
		cout<<"CHP MT for witness path"<<runT<<endl;
		writeShortCutorder(indexfile);
		writeCHPIncrease(indexfile);

		ReadShortCut(indexfile);
		readCHPIncrease(indexfile);
	}
	else
	{
		ReadShortCut(indexfile);
		readCHPIncrease(indexfile);
	}

	return 0;
}

//witness path retrieval
int Graph::CHcontractionorderMTWP(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult){
	benchmark::heap<2,int,int> Heap(nodenum);
	vector<int> vDistance(nodenum, INF);
	vector<int> vPre(nodenum,-1);
	int topNodeID, topDistance, neighborNodeID, neighborLength;
	int w;
	vDistance[ID1]=0;
	Heap.update(ID1,0);
	map<int,int> mWDistance;
	vector<pair<int,int>> vpWDistance;
	map<int,int> mDistance;
	int maxWDistance=-1;

	for(auto ivp=vW.begin();ivp!=vW.end();ivp++){
		w=(*ivp).first;
		int d=(*ivp).second+dUV;
		mWDistance[w]=d;
		if(d>maxWDistance)
			maxWDistance=d;
		mDistance[w]=INF;
	}

	if(mDistance.empty()){
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
				vPre[neighborNodeID]=topNodeID;
				if(mWDistance.find(neighborNodeID) != mWDistance.end())
					mDistance[neighborNodeID] = d;
			}else if(vDistance[neighborNodeID]>d){
				vDistance[neighborNodeID] = d;
				Heap.update(neighborNodeID, d);
				vPre[neighborNodeID]=topNodeID;
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

	//mDistance is the Dijkstra's distance; mWDistance is the sum of two weight
	for(auto imDistance = mDistance.begin(); imDistance != mDistance.end(); imDistance++)
	{
		int w = (*imDistance).first;
		if((*imDistance).second > mWDistance[(*imDistance).first])
		{
			//int w = (*imDistance).first;
			int distance = mWDistance[(*imDistance).first];

			vvpResult[k].push_back(make_pair(w, distance));
		}
		else{
			//cout<<"wp start: "<<ID1<<" "<<w<<endl;
			//retrieve the witness path (sketch path)
			vector<int> path;
			path.push_back(w);
			int pre=vPre[w];
			/*if(ID1==6644 && w==8310){
				cout<<"////// "<<pre<<endl;
			}*/
			while(pre!=ID1){
				path.insert(path.begin(), pre);
				pre=vPre[pre];
				/*if(ID1==6644 && w==8310){
					cout<<pre<<endl;
				}*/
			}
			path.insert(path.begin(), pre);

			//to get the detailed path
			for(int i=0;i<path.size()-1;i++){
				int currentID = path[i];
				while(OutEdgesM[currentID].find(path[i+1])!=OutEdgesM[currentID].end()){//find the original adjacent node
					int MnodeID=OutEdgesM[currentID][path[i+1]];
					path.insert(path.begin()+i+1, MnodeID);//add the MnodeID after the current i-th node
					/*if(ID1==6644 && w==8310){
						cout<<MnodeID<<" ,currentID: "<<currentID<<" "<<path[i+1]<<endl;
					}*/
				}
			}//path finish retrieval

			//recompute and check the path length
			/*cout<<"path length computating ";
			//path length check
			int plength=0;
			for(int i=0;i<path.size()-1;i++){
				int e1=path[i];
				int e2=path[i+1];
				for(int k=0;k<Neighbor[e1].size();k++){
					if(Neighbor[e1][k].first==e2){
						plength+=Neighbor[e1][k].second;
						break;
					}
				}
			}
			cout<<plength<<" , path calculated length "<<(*imDistance).second<<endl;*/

			//store the path information
			int a,b;
			for(int k=0;k<path.size()-1;k++){
				if(path[k]<path[k+1]){
					a=path[k];
					b=path[k+1];
				}else{
					a=path[k+1];
					b=path[k];
				}
				int edgeID=EdgeRe[make_pair(a,b)];
				vSm[edgeID]->wait();
				tri TRI;
				TRI.u=ID1;
				TRI.v=ID2;
				TRI.w=w;
				EdgeOnPath[edgeID].push_back(TRI);
				vSm[edgeID]->notify();
			}

			PathInfor[ID1].insert(make_pair(make_pair(ID2,w),(*imDistance).second));
			//cout<<"wp finish! "<<path.size()<<endl;

		}
	}

	return 0;
}

void Graph::CHPconsorderMTWP(){
	vector<pair<int,int>> nodeinfor;
	vvpShortCut.assign(nodenum, nodeinfor);
	vector<bool> vbVisited(nodenum, false);
	map<int, int> mPosition;
	mPosition[0]=nodenum-1;
	bool bUpdated;
	vector<pair<int,int>> vU,vW;
	map<int,int> mapinfor;
	OutEdgesM.assign(nodenum, mapinfor);
	/*unordered_map<pair<int,int>, int, hash_pair> pathinfvec;
	PathInfor.assign(nodenum, pathinfvec);
	vector<tri> vec;
	EdgeOnPath.assign(Edge.size(),vec);*/

	//cout<<"Building static CH "<<endl;

	int v;// current contracting vertex
	for(int i=0;i<vNodeOrder.size()-1;i++){
		//if(i%1000==0)
			//cout<<"i=/////////////////////// "<<i<<endl;
		v=vNodeOrder[i];
		if(v!=-1){
			vU=vvNode[v];
			//cout<<"contracting node "<<v<<endl;

			//delete the redundant ID1
			map<int,int> ID1s;
			for(int l=0;l<vU.size();l++){
				int ID1=vU[l].first;
				int wei1=vU[l].second;
				if(ID1s.find(ID1)==ID1s.end()){
					ID1s.insert(vU[l]);
				}else{
					if(ID1s[ID1]>wei1)
						ID1s[ID1]=wei1;
				}
			}

			vector<pair<int,int>> vU1;
			for(auto it=ID1s.begin();it!=ID1s.end();it++){
				vU1.push_back(*it);
			}

			//multiple thread
			vector<vector<pair<int,int>>> vvpResult;
			vector<pair<int,int>> vec;
			vec.clear();
			vvpResult.assign(vU1.size(), vec);
			//filter the neighbors no need for contraction
			vector<vector<pair<int,int>>> vUeach;
			vUeach.assign(vU1.size(), vec);
			int w;
			for(int l=0;l<vU1.size();l++){
				int ID1=vU1[l].first;
				for(auto ivp=vU1.begin();ivp!=vU1.end();ivp++){
					w=(*ivp).first;
					if(NodeOrder[w]>NodeOrder[ID1]){//to get rid of redundant computation
						if(vbVisited[w])
							continue;
						vUeach[l].push_back(*ivp);
						if(w<ID1)
							SupportNodes[w][ID1].push_back(v);
						else
							SupportNodes[ID1][w].push_back(v);
					}
				}
			}

			boost::thread_group threadf;
			for(int k=0;k<vU1.size();k++){
				if(vUeach[k].size()>0)
					threadf.add_thread(new boost::thread(&Graph::CHcontractionorderMTWP, this, k, vU1[k].first, v, boost::ref(vbVisited), vU1[k].second, boost::ref(vUeach[k]), boost::ref(vvpResult)));
			}
			threadf.join_all();
			//merge to new shortcut and write it back
			for(int k=0;k<vU1.size();k++)
			{
				int	ID1 = vU1[k].first;
				for(int j=0;j<vvpResult[k].size();j++){
					int w=vvpResult[k][j].first;
					int distance=vvpResult[k][j].second;
					vvpShortCut[ID1].push_back(make_pair(w, distance));
					vvNode[ID1].push_back(make_pair(w, distance));
					vvNode[w].push_back(make_pair(ID1, distance));
					//store the contracted node for path retrieval
					OutEdgesM[ID1].insert(make_pair(w,v));
					OutEdgesM[w].insert(make_pair(ID1,v));
				}
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

}

void Graph::insertEorder(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
		//DD[u]++;
		//DD2[u]++;
	}
	else{
		if(E[u][v].first>w)
			E[u][v]=make_pair(w,1);
		else if(E[u][v].first==w)
			E[u][v].second+=1;
	}

	if(E[v].find(u)==E[v].end()){
		E[v].insert(make_pair(u,make_pair(w,1)));
		//DD[v]++;
		//DD2[v]++;
	}
	else{
		if(E[v][u].first>w)
			E[v][u]=make_pair(w,1);
		else if(E[v][u].first==w)
			E[v][u].second+=1;
	}
}

void Graph::deleteEorder(int u,int v){
	if(E[u].find(v)!=E[u].end()){
		E[u].erase(E[u].find(v));
		//DD[u]--;
	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
		//DD[v]--;
	}
}



void Graph::deleteE(int u,int v){
	if(E[u].find(v)!=E[u].end()){
		E[u].erase(E[u].find(v));
		DD[u]--;
	}

	if(E[v].find(u)!=E[v].end()){
		E[v].erase(E[v].find(u));
		DD[v]--;
	}
}

vector<int> NodeOrders;
struct OrderComp{
	int x;
	int y;//order(x)<order(y)
	OrderComp(int _x, int _y){
		x=_x; y=_y;
	}
	bool operator< (const OrderComp& d) const{
		if(x==d.x && y==d.y){//avoid the redundant
			return false;
		}else{
			if(x!=d.x)
				return NodeOrders[x]<NodeOrders[d.x];
			if(y!=d.y)
				return NodeOrders[y]<NodeOrders[d.y];
		}
	}
};

int Graph::NewSCweight(int s, int t){//new weight of shortcut
	int wt=INF;
	for(int i=0;i<Neighbor[s].size();i++){
		if(Neighbor[s][i].first==t){
			wt=Neighbor[s][i].second;//the weight value in the original graph
			break;
		}
	}
	int ssw,wtt,wid;
	vector<int> Wnodes; //Wnodes.clear();
	if(s<t)
		Wnodes=SupportNodes[s][t];
	else
		Wnodes=SupportNodes[t][s];
	for(int i=0;i<Wnodes.size();i++){
		wid=Wnodes[i];
		for(int j=0;j<AdjaShort[wid].size();j++){
			if(AdjaShort[wid][j].first==s){
				ssw=AdjaShort[wid][j].second;
			}
			if(AdjaShort[wid][j].first==t){
				wtt=AdjaShort[wid][j].second;
			}
		}
		if(ssw+wtt<wt){
			wt=ssw+wtt;
		}
	}
	return wt;
}

//the increase case of pruned CH
void Graph::CHPinc(int a, int b, int oldW, int newW){
	//modify the original graph information
	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int i=0;i<Neighbor[b].size();i++){
		if(Neighbor[b][i].first==a){
			Neighbor[b][i].second=newW;
			break;
		}
	}

	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC; //OC.clear();
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	//refresh the weight of existing shortcuts
	if(NodeOrder[a]<NodeOrder[b]){
		for(int i=0;i<AdjaShort[a].size();i++){
			if(AdjaShort[a][i].first==b){
				if(AdjaShort[a][i].second==oldW){
					OrderComp oc={a,b};
					OC.insert(oc);
					OCdis[make_pair(a,b)]=oldW;
				}
				break;
			}
		}
	}else{
		for(int i=0;i<AdjaShort[b].size();i++){
			if(AdjaShort[b][i].first==a){
				if(AdjaShort[b][i].second==oldW){
					OrderComp oc={b,a};
					OC.insert(oc);
					OCdis[make_pair(b,a)]=oldW;
				}
				break;
			}
		}
	}
	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];//distance of s--->t before change
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		//HigherIn.clear(); LowerIn.clear();
		//the shortcuts infected by s-->t
		for(int i=0;i<AdjaShort[s].size();i++){
			inID=AdjaShort[s][i].first;
			inW=AdjaShort[s][i].second;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<AdjaShort[t].size();i++){
			inID=AdjaShort[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(AdjaShort[t][i].second==wt+inW){
					OrderComp oc={t,inID};
					OC.insert(oc);
					OCdis[make_pair(t,inID)]=wt+inW;
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<AdjaShort[inID].size();j++){
				if(AdjaShort[inID][j].first==t){
					if(AdjaShort[inID][j].second==inW+wt){
						OrderComp oc={inID,t};
						OC.insert(oc);
						OCdis[make_pair(inID,t)]=wt+inW;
					}
					break;
				}
			}
		}

		//get the new weight value of s-->t
		wt=NewSCweight(s,t);

		for(int i=0;i<AdjaShort[s].size();i++){
			if(AdjaShort[s][i].first==t){
				AdjaShort[s][i].second=wt;
				break;
			}
		}
		for(int i=0;i<AdjaShortR[t].size();i++){
			if(AdjaShortR[t][i].first==s){
				AdjaShortR[t][i].second=wt;
				break;
			}
		}
	}

	//cout<<"existing shortcut finish updating!"<<endl;

	//check which shortcuts needs to be added through witness paths
	//1. find the witness path with (a,b) on it
	int lid, hid;
	if(a<b){
		lid=a;
		hid=b;
	}else{
		lid=b;
		hid=a;
	}
	int edgeID=EdgeRe[make_pair(lid,hid)];

	int pID, u, w, length, newlength;
	set<OrderComp> OCnew;
	map<pair<int,int>,int> OCdisnew;

	//cout<<"the number of witness path passing through "<<EdgeOnPath[EdgeRe[make_pair(lid,hid)]].size()<<endl;

	//find those witness paths passing through (a,b)
	for(int i=0;i<EdgeOnPath[edgeID].size();i++){
		tri TRI=EdgeOnPath[edgeID][i];
		if(InvalidWP[TRI.u].find(make_pair(TRI.v, TRI.w))==InvalidWP[TRI.u].end()){
			length=PathInfor[TRI.u][make_pair(TRI.v,TRI.w)];
			u=TRI.u;
			w=TRI.w;
			newlength=length+newW-oldW;
			//cout<<i<<" "<<u<<" "<<TRI.v<<" "<<w<<" "<<newlength<<endl;

			int wt=NewSCweight(u,w);
			if(wt<newlength){
				//cout<<"needs inserting new shortcut! "<<u<<" "<<w<<" "<<wt<<endl;
				OrderComp oc={u,w};
				OCnew.insert(oc);
				OCdisnew[make_pair(u,w)]=wt;
				InvalidWP[u].insert(make_pair(TRI.v,w));
			}else
				PathInfor[TRI.u][make_pair(TRI.v,TRI.w)]=newlength;//update the witness path length
		}
	}

	//insert the new shortcuts caused by the invalidation of witness paths
	while(!OCnew.empty()){
		int s=(*OCnew.begin()).x; int t=(*OCnew.begin()).y;
		OCnew.erase(OCnew.begin());
		int wt;
		wt=OCdisnew[make_pair(s,t)];

		//1-N Dijstra's search from t
		//maintain the SupportNodes at the same time
		//insert the newly added shortcuts into OCnew
		vector<pair<int,int>> addedShortcut;
		CHcontractInc(t,s,wt,AdjaShort[s],addedShortcut);
		int id, dis;
		//cout<<"addedShortcut size "<<addedShortcut.size()<<endl;
		for(int k=0;k<addedShortcut.size();k++){
			id=addedShortcut[k].first;
			dis=addedShortcut[k].second;
			if(NodeOrder[t]<NodeOrder[id]){
				OrderComp oc={t,id};
				OCnew.insert(oc);
				OCdisnew[make_pair(t,id)]=dis;
			}else{
				OrderComp oc={id,t};
				OCnew.insert(oc);
				OCdisnew[make_pair(id,t)]=dis;
			}
		}

		for(int x=0;x<AdjaShort[s].size();x++){
			int ID=AdjaShort[s][x].first;
			if(ID<t)
				SupportNodes[ID][t].push_back(s);
			else
				SupportNodes[t][ID].push_back(s);
		}
		AdjaShort[s].push_back(make_pair(t,wt));
		AdjaShortR[t].push_back(make_pair(s,wt));
	}

}

int Graph::CHcontractInc(int ID1, int ID2, int dUV, vector<pair<int, int> >& vW, vector<pair<int,int>>& addedShortcut){
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
		if(w!=ID1){
			int d=(*ivp).second+dUV;
			mWDistance[w]=d;
			if(d>maxWDistance)
				maxWDistance=d;
			mDistance[w]=INF;
		}
	}
	if(mDistance.empty()){
		return 0;
	}
	for(auto imWDistance=mWDistance.begin();imWDistance!=mWDistance.end();imWDistance++)
		vpWDistance.push_back(make_pair((*imWDistance).first, (*imWDistance).second));

	int dThreshold=maxWDistance;

	while(!Heap.empty()){
		Heap.extract_min(topNodeID, topDistance);
		//if(vbVisited[topNodeID])
			//continue;
		if(topDistance>dThreshold)
			break;
		//to higher neighbors
		for(auto ivp=AdjaShort[topNodeID].begin();ivp!=AdjaShort[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(NodeOrder[neighborNodeID]<=NodeOrder[ID2])
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
		//to lower neighbors
		for(auto ivp=AdjaShortR[topNodeID].begin();ivp!=AdjaShortR[topNodeID].end();ivp++){
			neighborNodeID = (*ivp).first;
			neighborLength = (*ivp).second;
			if(NodeOrder[neighborNodeID]<=NodeOrder[ID2])
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
		if((*imDistance).second > mWDistance[(*imDistance).first])
		{
			int w = (*imDistance).first;
			int distance = mWDistance[(*imDistance).first];

			addedShortcut.push_back(make_pair(w,distance));
		}
	}
	return 0;
}

//the decrease case of CH with pruning
void Graph::CHPdec(int a, int b, int oldW, int newW){
	//modify the original graph information
	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int i=0;i<Neighbor[b].size();i++){
		if(Neighbor[b][i].first==a){
			Neighbor[b][i].second=newW;
			break;
		}
	}
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the fresh distance and avoid search in the adjacent list
	//OC.clear(); OCdis.clear();

	if(NodeOrder[a]<NodeOrder[b]){
		for(int i=0;i<AdjaShort[a].size();i++){
			if(AdjaShort[a][i].first==b){
				if(AdjaShort[a][i].second>newW){
					AdjaShort[a][i].second=newW;

					OCdis[make_pair(a,b)]=newW;
					OC.insert(OrderComp(a,b));
				}
				break;
			}
		}
	}else{
		for(int i=0;i<AdjaShort[b].size();i++){
			if(AdjaShort[b][i].first==a){
				if(AdjaShort[b][i].second>newW){
					AdjaShort[b][i].second=newW;

					OCdis[make_pair(b,a)]=newW;
					OC.insert(OrderComp(b,a));
				}
				break;
			}
		}
	}

	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; //InM2t.clear();
		vector<pair<int,int>> InMLower; //InMLower.clear();
		for(int i=0;i<AdjaShort[s].size();i++){
			if(NodeOrder[AdjaShort[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(AdjaShort[s][i].first,AdjaShort[s][i].second));
			else if(NodeOrder[AdjaShort[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(AdjaShort[s][i].first,AdjaShort[s][i].second));
		}
		int inID,inW,inWt;
		for(int i=0;i<AdjaShort[t].size();i++){
			inID=AdjaShort[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=AdjaShort[t][i].second;
				if(inWt>inW+wt){
					AdjaShort[t][i].second=inW+wt;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<AdjaShort[inID].size();j++){
				if(AdjaShort[inID][j].first==t){
					inWt=AdjaShort[inID][j].second;
					if(inWt>inW+wt){
						AdjaShort[inID][j].second=inW+wt;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}
					break;
				}
			}
		}
	}//finish change index
}

void Graph::CHdecStr(int a,int b, int oldW, int newW){
	//modify the original graph information
	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int i=0;i<Neighbor[b].size();i++){
		if(Neighbor[b][i].first==a){
			Neighbor[b][i].second=newW;
			break;
		}
	}
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the fresh distance and avoid search in the adjacent list
	//OC.clear(); OCdis.clear();

	if(NodeOrder[a]<NodeOrder[b]){
		for(int i=0;i<NeighborCon[a].size();i++){
			if(NeighborCon[a][i].first==b){
				if(NeighborCon[a][i].second.first>newW){
					NeighborCon[a][i].second.first=newW;
					NeighborCon[a][i].second.second=1;

					OCdis[make_pair(a,b)]=newW;
					OC.insert(OrderComp(a,b));
				}else if(NeighborCon[a][i].second.first==newW)
					NeighborCon[a][i].second.second+=1;
				break;
			}
		}
	}else{
		for(int i=0;i<NeighborCon[b].size();i++){
			if(NeighborCon[b][i].first==a){
				if(NeighborCon[b][i].second.first>newW){
					NeighborCon[b][i].second.first=newW;
					NeighborCon[b][i].second.second=1;

					OCdis[make_pair(b,a)]=newW;
					OC.insert(OrderComp(b,a));
				}else if(NeighborCon[b][i].second.first==newW)
					NeighborCon[b][i].second.second+=1;
				break;
			}
		}
	}


	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; //InM2t.clear();
		vector<pair<int,int>> InMLower; //InMLower.clear();
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
			else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
		}
		int inID,inW,inWt;
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=NeighborCon[t][i].second.first;
				if(inWt>inW+wt){
					NeighborCon[t][i].second.first=inW+wt;
					NeighborCon[t][i].second.second=1;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}else if(inWt==inW+wt){
					NeighborCon[t][i].second.second+=1;
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					inWt=NeighborCon[inID][j].second.first;
					if(inWt>inW+wt){
						NeighborCon[inID][j].second.first=inW+wt;
						NeighborCon[inID][j].second.second=1;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}else if(inWt==inW+wt)
						NeighborCon[inID][j].second.second+=1;
					break;
				}
			}
		}
	}//finish change index
}

void Graph::CHincStrMT(int a,int b, int oldW, int newW){
	//modify the original graph information
	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int i=0;i<Neighbor[b].size();i++){
		if(Neighbor[b][i].first==a){
			Neighbor[b][i].second=newW;
			break;
		}
	}

	//NodeOrders.clear();
		NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
		set<OrderComp> OC; //OC.clear();
		map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
		//OCdis.clear();

		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}


		while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());

		wt=OCdis[make_pair(s,t)];//distance of s--->t before change
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		//HigherIn.clear(); LowerIn.clear();
		//the shortcuts infected by s-->t
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		//get the new weight value of s-->t
		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;//the weight value in the original graph
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; //Wnodes.clear();
		if(s<t){
			//Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
			Wnodes=SCconNodesMT[s][t];
		}else{
			//Wnodes=SCconNodes[make_pair(t,s)];
			Wnodes=SCconNodesMT[t][s];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}

		//refresh the weight value of s--t in the index
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){
				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}
}

void Graph::CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	//maintain the index caused by the weight change
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC;
	map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
	//OC.clear(); OCdis.clear();

	int a,b,newW;//the weight of (a,b) decrease to newW
	for(int k=0;k<wBatch.size();k++){
		a=wBatch[k].first.first;
		b=wBatch[k].first.second;
		newW=wBatch[k].second.second;

		//modify the information in original graph
		for(int i=0;i<Neighbor[a].size();i++){
			if(Neighbor[a][i].first==b){
				Neighbor[a][i].second=newW;
				break;
			}
		}
		for(int i=0;i<Neighbor[b].size();i++){
			if(Neighbor[b][i].first==a){
				Neighbor[b][i].second=newW;
				break;
			}
		}

		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first>newW){
						//cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
						NeighborCon[a][i].second.first=newW;
						NeighborCon[a][i].second.second=1;

						OCdis[make_pair(a,b)]=newW;
						OC.insert(OrderComp(a,b));
					}else if(NeighborCon[a][i].second.first==newW)
						NeighborCon[a][i].second.second+=1;
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first>newW){
						NeighborCon[b][i].second.first=newW;
						NeighborCon[b][i].second.second=1;

						OCdis[make_pair(b,a)]=newW;
						OC.insert(OrderComp(b,a));
					}else if(NeighborCon[b][i].second.first==newW)
						NeighborCon[b][i].second.second+=1;
					break;
				}
			}
		}
	}


	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());
		wt=OCdis[make_pair(s,t)];
		map<int,int> InM2t; //InM2t.clear();
		vector<pair<int,int>> InMLower; //InMLower.clear();
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
				InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
			else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
				InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
		}
		int inID,inW,inWt;
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(InM2t.find(inID)!=InM2t.end()){
				inW=InM2t[inID];
				inWt=NeighborCon[t][i].second.first;
				if(inWt>inW+wt){
					NeighborCon[t][i].second.first=inW+wt;
					NeighborCon[t][i].second.second=1;
					OCdis[make_pair(t,inID)]=inW+wt;
					OrderComp oc={t,inID};
					OC.insert(oc);
				}else if(inWt==inW+wt){
					NeighborCon[t][i].second.second+=1;
				}
			}
		}

		for(int i=0;i<InMLower.size();i++){
			inID=InMLower[i].first; inW=InMLower[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					inWt=NeighborCon[inID][j].second.first;
					if(inWt>inW+wt){
						NeighborCon[inID][j].second.first=inW+wt;
						NeighborCon[inID][j].second.second=1;

						OCdis[make_pair(inID,t)]=inW+wt;
						OrderComp oc={inID,t};
						OC.insert(oc);
					}else if(inWt==inW+wt)
						NeighborCon[inID][j].second.second+=1;
					break;
				}
			}
		}
	}//finish change index
}

void Graph::CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	//NodeOrders.clear();
	NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
	set<OrderComp> OC; //OC.clear();
	map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	for(int wb=0;wb<wBatch.size();wb++){
		int a=wBatch[wb].first.first;
		int b=wBatch[wb].first.second;
		int oldW=wBatch[wb].second.first;
		int newW=wBatch[wb].second.second;

		//modify the original graph information
		for(int i=0;i<Neighbor[a].size();i++){
			if(Neighbor[a][i].first==b){
				Neighbor[a][i].second=newW;
				break;
			}
		}
		for(int i=0;i<Neighbor[b].size();i++){
			if(Neighbor[b][i].first==a){
				Neighbor[b][i].second=newW;
				break;
			}
		}


		if(NodeOrder[a]<NodeOrder[b]){
			for(int i=0;i<NeighborCon[a].size();i++){
				if(NeighborCon[a][i].first==b){
					if(NeighborCon[a][i].second.first==oldW){
						NeighborCon[a][i].second.second-=1;
						if(NeighborCon[a][i].second.second<1){
							OrderComp oc={a,b};
							OC.insert(oc);
							OCdis[make_pair(a,b)]=oldW;
						}
					}
					break;
				}
			}
		}else{
			for(int i=0;i<NeighborCon[b].size();i++){
				if(NeighborCon[b][i].first==a){
					if(NeighborCon[b][i].second.first==oldW){
						NeighborCon[b][i].second.second-=1;
						if(NeighborCon[b][i].second.second<1){
							OrderComp oc={b,a};
							OC.insert(oc);
							OCdis[make_pair(b,a)]=oldW;
						}
					}
					break;
				}
			}
		}
	}

	while(!OC.empty()){
		int s=(*OC.begin()).x; int t=(*OC.begin()).y;
		int wt;
		OC.erase(OC.begin());
		wt=OCdis[make_pair(s,t)];//distance of s--->t before change
		int inID,inW;
		map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
		//HigherIn.clear(); LowerIn.clear();
		//the shortcuts infected by s-->t
		for(int i=0;i<NeighborCon[s].size();i++){
			inID=NeighborCon[s][i].first;
			inW=NeighborCon[s][i].second.first;
			if(NodeOrder[inID]<NodeOrder[t]){
				LowerIn.push_back(make_pair(inID,inW));
			}else if(NodeOrder[inID]>NodeOrder[t]){
				HigherIn.insert(make_pair(inID,inW));
			}
		}
		for(int i=0;i<NeighborCon[t].size();i++){
			inID=NeighborCon[t][i].first;
			if(HigherIn.find(inID)!=HigherIn.end()){
				inW=HigherIn[inID];
				if(NeighborCon[t][i].second.first==wt+inW){
					NeighborCon[t][i].second.second-=1;
					if(NeighborCon[t][i].second.second<1){
						OrderComp oc={t,inID};
						OC.insert(oc);
						OCdis[make_pair(t,inID)]=wt+inW;
					}
				}
			}
		}
		for(int i=0;i<LowerIn.size();i++){
			inID=LowerIn[i].first; inW=LowerIn[i].second;
			for(int j=0;j<NeighborCon[inID].size();j++){
				if(NeighborCon[inID][j].first==t){
					if(NeighborCon[inID][j].second.first==inW+wt){
						NeighborCon[inID][j].second.second-=1;
						if(NeighborCon[inID][j].second.second<1){
							OrderComp oc={inID,t};
							OC.insert(oc);
							OCdis[make_pair(inID,t)]=wt+inW;
						}
					}
					break;
				}
			}
		}

		//get the new weight value of s-->t
		wt=INF; int countwt=0;
		for(int i=0;i<Neighbor[s].size();i++){
			if(Neighbor[s][i].first==t){
				wt=Neighbor[s][i].second;//the weight value in the original graph
				countwt=1;
				break;
			}
		}
		int ssw,wtt,wid;
		vector<int> Wnodes; //Wnodes.clear();
		if(s<t){
			//Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
			Wnodes=SCconNodesMT[s][t];
		}else{
			//Wnodes=SCconNodes[make_pair(t,s)];
			Wnodes=SCconNodesMT[s][t];
		}

		for(int i=0;i<Wnodes.size();i++){
			wid=Wnodes[i];
			for(int j=0;j<NeighborCon[wid].size();j++){
				if(NeighborCon[wid][j].first==s){
					ssw=NeighborCon[wid][j].second.first;
				}
				if(NeighborCon[wid][j].first==t){
					wtt=NeighborCon[wid][j].second.first;
				}
			}

			if(ssw+wtt<wt){
				wt=ssw+wtt;
				countwt=1;
			}else if(ssw+wtt==wt){
				countwt+=1;
			}
		}

		//refresh the weight value of s--t in the index
		for(int i=0;i<NeighborCon[s].size();i++){
			if(NeighborCon[s][i].first==t){
				NeighborCon[s][i].second.first=wt;
				NeighborCon[s][i].second.second=countwt;
				break;
			}
		}
	}
}
