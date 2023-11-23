/*
 * Labelcon.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::H2HconOrderMT(string orderfile){
    ifstream IF(orderfile);
    if(!IF){
        cout<<"Cannot open Map "<<orderfile<<endl;
        CHconsMTOrderGenerate(orderfile);//generating vertex ordering by minimum degree elimination
        exit(0);
    }else{
        IF.close();
        CHconsorderMT(orderfile);
    }
	makeTree();
	makeIndex();
}

void Graph::makeTree(){
	vector<int> vecemp; //vecemp.clear();
	VidtoTNid.assign(nodenum,vecemp);

	rank.assign(nodenum,0);
	//Tree.clear();
	int len=vNodeOrder.size()-1;
	heightMax=0;

	Node rootn;
	int x=vNodeOrder[len];
	//cout<<"len "<<len<<" , ID "<<x<<endl;
	while(x==-1){//to skip those vertices whose ID is -1
		len--;
		x=vNodeOrder[len];
		//cout<<"len "<<len<<" , ID "<<x<<endl;
	}
	rootn.vert=NeighborCon[x];
	rootn.uniqueVertex=x;
	rootn.pa=-1;
	rootn.height=1;
	rank[x]=0;
	Tree.push_back(rootn);
	len--;

	int nn;
	for(;len>=0;len--){
		int x=vNodeOrder[len];
		Node nod;
		nod.vert=NeighborCon[x];
		nod.uniqueVertex=x;
		int pa=match(x,NeighborCon[x]);
		Tree[pa].ch.push_back(Tree.size());
		nod.pa=pa;
		nod.height=Tree[pa].height+1;

		nod.hdepth=Tree[pa].height+1;
		for(int i=0;i<NeighborCon[x].size();i++){
			nn=NeighborCon[x][i].first;
			VidtoTNid[nn].push_back(Tree.size());
			if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
				Tree[rank[nn]].hdepth=Tree[pa].height+1;
		}
		if(nod.height>heightMax) heightMax=nod.height;
		rank[x]=Tree.size();
		Tree.push_back(nod);
		//cout<<"len "<<len<<" , ID "<<x<<endl;
	}
}

void Graph::makeIndex(){
	makeRMQ();

	//initialize
	vector<int> list; //list.clear();
	list.push_back(Tree[0].uniqueVertex);
	Tree[0].pos.clear();
	Tree[0].pos.push_back(0);

	for(int i=0;i<Tree[0].ch.size();i++){
		makeIndexDFS(Tree[0].ch[i],list);
	}

}

int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
	int nearest=vert[0].first;
	for(int i=1;i<vert.size();i++){
		if(rank[vert[i].first]>rank[nearest])
			nearest=vert[i].first;
	}
	int p=rank[nearest];
	return p;
}

void Graph::makeRMQDFS(int p, int height){
	toRMQ[p] = EulerSeq.size();
	EulerSeq.push_back(p);
	for (int i = 0; i < Tree[p].ch.size(); i++){
		makeRMQDFS(Tree[p].ch[i], height + 1);
		EulerSeq.push_back(p);
	}
}

void Graph::makeRMQ(){
	//EulerSeq.clear();
	toRMQ.assign(nodenum,0);
	//RMQIndex.clear();
	makeRMQDFS(0, 1);
	RMQIndex.push_back(EulerSeq);

	int m = EulerSeq.size();
	for (int i = 2, k = 1; i < m; i = i * 2, k++){
		vector<int> tmp;
		//tmp.clear();
		tmp.assign(m,0);
		for (int j = 0; j < m - i; j++){
			int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
			if (Tree[x].height < Tree[y].height)
				tmp[j] = x;
			else tmp[j] = y;
		}
		RMQIndex.push_back(tmp);
	}
}

int Graph::LCAQuery(int _p, int _q){
	int p = toRMQ[_p], q = toRMQ[_q];
	if (p > q){
		int x = p;
		p = q;
		q = x;
	}
	int len = q - p + 1;
	int i = 1, k = 0;
	while (i * 2 < len){
		i *= 2;
		k++;
	}
	q = q - i + 1;
	if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
		return RMQIndex[k][p];
	else return RMQIndex[k][q];
}

void Graph::makeIndexDFS(int p, vector<int>& list){
	//initialize
	int NeiNum=Tree[p].vert.size();
	Tree[p].pos.assign(NeiNum+1,0);
	Tree[p].dis.assign(list.size(),INF);
	Tree[p].cnt.assign(list.size(),0);
	Tree[p].FN.assign(list.size(),true);
    Tree[p].vAncestor=list;
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
				}else if(Tree[p].dis[j]==z+disvb){
					Tree[p].cnt[j]+=1;
				}
			}
		}
	}

	//nested loop
	list.push_back(Tree[p].uniqueVertex);
	for(int i=0;i<Tree[p].ch.size();i++){
		makeIndexDFS(Tree[p].ch[i],list);
	}
	list.pop_back();
}

vector<int> NodeOrderss;
struct OrderCompp{//prior to reture the vertex with smaller order
	int x;
	OrderCompp(int _x){
		x=_x;
	}
	bool operator< (const OrderCompp& d) const{
		if(x==d.x){//avoid the redundant
			return false;
		}else{
			if(x!=d.x)
				return NodeOrderss[x]<NodeOrderss[d.x];
		}
	}
};

void Graph::H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	map<int,int> checkedDis;

	for(int i=0;i<Tree.size();i++){
		Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
	}

	//NodeOrderss.clear();
	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<int>> SCre; //SCre.clear();
	set<int> ss; //ss.clear();
	SCre.assign(nodenum,ss);//{vertexID, set<int>}
	set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

	set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

	int a,b,oldW,newW,lid,hid;
	for(int k=0;k<wBatch.size();k++){
		a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
		if(NodeOrder[a]<NodeOrder[b]){
			lid=a;hid=b;
		}else{
			lid=b;hid=a;
		}

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

		for(int i=0;i<Tree[rank[lid]].vert.size();i++){
			if(Tree[rank[lid]].vert[i].first==hid){
				if(Tree[rank[lid]].vert[i].second.first>newW){
					Tree[rank[lid]].vert[i].second.first=newW;
					Tree[rank[lid]].vert[i].second.second=1;
					SCre[lid].insert(hid);
					OC.insert(OrderCompp(lid));
				}else if(Tree[rank[lid]].vert[i].second.first==newW){
					Tree[rank[lid]].vert[i].second.second+=1;
				}
				break;
			}
		}

	}

	vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
	vector<int> ProBeginVertexSetNew;
	int ProBeginVertexID;
	int ProID;
	//processing the stars
	while(!OC.empty()){
		ProID=(*OC.begin()).x;
		OC.erase(OC.begin());
		vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
		bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
		for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
			int Cid=*it; int Cw;
			int cidH=Tree[rank[Cid]].height-1;

			map<int,int> Hnei; //Hnei.clear();
			vector<pair<int,int>> Lnei; //Lnei.clear();
			for(int j=0;j<Vert.size();j++){
				if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
					Hnei[Vert[j].first]=Vert[j].second.first;
				}else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
					Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
				}else{
					Cw=Vert[j].second.first;
				}
			}

			if(Tree[rank[ProID]].dis[cidH]>Cw){
				Tree[rank[ProID]].dis[cidH]=Cw;
				Tree[rank[ProID]].FN[cidH]=true;
				ProIDdisCha=true;
				Tree[rank[ProID]].DisRe.insert(Cid);
			}else if(Tree[rank[ProID]].dis[cidH]==Cw){
				Tree[rank[ProID]].FN[cidH]=true;
			}

			int hid,hidHeight,lid,lidHeight,wsum;
			for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
				hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
				if(Hnei.find(hid)!=Hnei.end()){
					wsum=Cw+Hnei[hid];
					if(wsum<Tree[rank[Cid]].vert[j].second.first){
						Tree[rank[Cid]].vert[j].second.first=wsum;
						Tree[rank[Cid]].vert[j].second.second=1;
						SCre[Cid].insert(hid);
						OC.insert(OrderCompp(Cid));
					}else if(wsum==Tree[rank[Cid]].vert[j].second.first){
						Tree[rank[Cid]].vert[j].second.second+=1;
					}

				}
			}
			for(int j=0;j<Lnei.size();j++){
				lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
				for(int k=0;k<Tree[rank[lid]].vert.size();k++){
					if(Tree[rank[lid]].vert[k].first==Cid){
						wsum=Cw+Lnei[j].second;
						if(Tree[rank[lid]].vert[k].second.first>wsum){
							Tree[rank[lid]].vert[k].second.first=wsum;
							Tree[rank[lid]].vert[k].second.second=1;
							SCre[lid].insert(Cid);
							OC.insert(OrderCompp(lid));
						}else if(Tree[rank[lid]].vert[k].second.first==wsum){
							Tree[rank[lid]].vert[k].second.second+=1;
						}

						break;
					}
				}
			}
		}

		if(ProIDdisCha){//if the distance labeling is dectected changed
			vertexIDChL.insert(ProID);
			ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
			ProBeginVertexSetNew.push_back(ProID);
			int rnew=rank[ProID],r;
			for(int i=0;i<ProBeginVertexSet.size();i++){
				r=rank[ProBeginVertexSet[i]];
				if(LCAQuery(rnew,r)!=rnew){
					ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
				}
			}
			ProBeginVertexSet=ProBeginVertexSetNew;
		}
	}

	//cout<<"Finish bottom-up refresh"<<endl;
	for(int i=0;i<ProBeginVertexSet.size();i++){
		ProBeginVertexID=ProBeginVertexSet[i];
		vector<int> linee; //linee.clear();
		linee.reserve(heightMax);
		int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
		while(Tree[rank[pachidd]].height>1){
			linee.insert(linee.begin(),pachidd);
			pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
		}
		linee.insert(linee.begin(),pachidd);
		EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis);
	}
	//return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
	bool ProIDdisCha=false;

	if(Tree[child].DisRe.size()!=0){
		for(int k=0;k<Tree[child].vert.size();k++){
			int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
			if(Tree[child].FN[bH]){
				if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
					for(int i=0;i<bH;i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
							Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
							Tree[child].FN[i]=false;
							ProIDdisCha=true;
						}
					}
					for(int i=bH+1;i<line.size();i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
							Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
							Tree[child].FN[i]=false;
							ProIDdisCha=true;
						}
					}

				}else{//partial ancestor check

					if(vertexIDChL.find(b)!=vertexIDChL.end()){
						for(int i=0;i<bH;i++){
							checkedDis.insert(make_pair(child,i));
							if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
								Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
								Tree[child].FN[i]=false;
								ProIDdisCha=true;
							}
						}
					}
					for(int i=bH+1;i<line.size();i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
							Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
							Tree[child].FN[i]=false;
							ProIDdisCha=true;
						}
					}

				}
			}
		}
	}else{
		for(int k=0;k<Tree[child].vert.size();k++){
			int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
			if(Tree[child].FN[bH]){
				if(vertexIDChL.find(b)!=vertexIDChL.end()){
					for(int i=0;i<bH;i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
							Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
							Tree[child].FN[i]=false;
							ProIDdisCha=true;
						}
					}
				}
				for(int i=bH+1;i<line.size();i++){
					checkedDis.insert(make_pair(child,i));
					if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
						Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
						Tree[child].FN[i]=false;
						ProIDdisCha=true;
					}
				}
			}
		}
	}

	if(ProIDdisCha){
		vertexIDChL.insert(Tree[child].uniqueVertex);
	}

	line.push_back(Tree[child].uniqueVertex);
	for(int i=0;i<Tree[child].ch.size();i++){
		EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis);
	}
	line.pop_back();

}

void Graph::H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	int checknum=0;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	//NodeOrderss.clear();
	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<int>> SCre; //SCre.clear();
	set<int> ss; //ss.clear();
	SCre.assign(nodenum,ss);//{vertexID, set<int>}
	set<OrderCompp> OC; OC.clear();//vertexID in decreasing node order

	for(int k=0;k<wBatch.size();k++){
		int a=wBatch[k].first.first;
		int b=wBatch[k].first.second;
		int oldW=wBatch[k].second.first;
		int newW=wBatch[k].second.second;

		if(oldW!=newW){
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

		int lid,hid;
		if(NodeOrder[a]<NodeOrder[b]){
			lid=a;hid=b;
		}else{
			lid=b;hid=a;
		}

		for(int i=0;i<Tree[rank[lid]].vert.size();i++){
			if(Tree[rank[lid]].vert[i].first==hid){
				if(Tree[rank[lid]].vert[i].second.first==oldW){
					Tree[rank[lid]].vert[i].second.second-=1;
					if(Tree[rank[lid]].vert[i].second.second<1){
						OCdis[make_pair(lid,hid)]=oldW;
						SCre[lid].insert(hid);
						OC.insert(OrderCompp(lid));
					}
				}
				break;
			}
		}
	}
	}

	vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
	vector<int> ProBeginVertexSetNew;
	bool influence;
	int ProID; vector<int> line;
	while(!OC.empty()){
		ProID=(*OC.begin()).x;
		OC.erase(OC.begin());
		vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
		influence=false;

		//each ProID corresponds to a line
		line.clear(); line.reserve(heightMax);
		int pachid=ProID;
		while(Tree[rank[pachid]].height>1){
			line.insert(line.begin(),pachid);
			pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
		}
		line.insert(line.begin(),pachid);

		for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
			int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
			int cidH=Tree[rank[Cid]].height-1;

			map<int,int> Hnei; //Hnei.clear();
			vector<pair<int,int>> Lnei; //Lnei.clear();
			for(int j=0;j<Vert.size();j++){
				if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
					Hnei[Vert[j].first]=Vert[j].second.first;
				}else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
					Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
				}
			}
			//check the affected shortcuts
			int hid,lid;
			for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
				hid=Tree[rank[Cid]].vert[j].first;
				if(Hnei.find(hid)!=Hnei.end()){
					if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
						Tree[rank[Cid]].vert[j].second.second-=1;
						if(Tree[rank[Cid]].vert[j].second.second<1){
							SCre[Cid].insert(hid);
							OC.insert(OrderCompp(Cid));
							OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
						}
					}
				}
			}
			for(int j=0;j<Lnei.size();j++){
				lid=Lnei[j].first;
				for(int k=0;k<Tree[rank[lid]].vert.size();k++){
					if(Tree[rank[lid]].vert[k].first==Cid){
						if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
							Tree[rank[lid]].vert[k].second.second-=1;
							if(Tree[rank[lid]].vert[k].second.second<1){
								SCre[lid].insert(Cid);
								OC.insert(OrderCompp(lid));
								OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
							}
						}
						break;
					}
				}
			}


			//before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
			if(Tree[rank[ProID]].FN[cidH]){
				influence=true;
				//higher than Cid
				for(int i=0;i<cidH;i++){
					if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
						Tree[rank[ProID]].cnt[i]-=1;
					}
				}

				//equal to Cid
				Tree[rank[ProID]].FN[cidH]=false;
				Tree[rank[ProID]].cnt[cidH]-=1;

				//lower than Cid
				for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
					if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
						Tree[rank[ProID]].cnt[i]-=1;
					}
				}
			}

			//get the new value of shortcut
		//	cout<<Cw<<" increase to ";
			Cw=INF; int countwt=0;

			for(int i=0;i<Neighbor[ProID].size();i++){
				if(Neighbor[ProID][i].first==Cid){
					Cw=Neighbor[ProID][i].second;//the weight value in the original graph
					countwt=1;
					break;
				}
			}

			int ssw,wtt,wid;
			vector<int> Wnodes;
			Wnodes.clear();

			if(ProID<Cid)
				Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
			else
				Wnodes=SCconNodesMT[Cid][ProID];
			if(Wnodes.size()>0){
				for(int i=0;i<Wnodes.size();i++){
					wid=Wnodes[i];
					for(int j=0;j<Tree[rank[wid]].vert.size();j++){
						if(Tree[rank[wid]].vert[j].first==ProID){
							ssw=Tree[rank[wid]].vert[j].second.first;
						}
						if(Tree[rank[wid]].vert[j].first==Cid){
							wtt=Tree[rank[wid]].vert[j].second.first;
						}
					}

					if(ssw+wtt<Cw){
						Cw=ssw+wtt;
						countwt=1;
					}else if(ssw+wtt==Cw){
						countwt+=1;
					}
				}
			}

			//cout<<Cw<<endl;
			//refresh the shortcut to the new value
			for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
				if(Tree[rank[ProID]].vert[i].first==Cid){
					Tree[rank[ProID]].vert[i].second.first=Cw;
					Tree[rank[ProID]].vert[i].second.second=countwt;
					break;
				}
			}
		}

		if(influence){
			ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
			ProBeginVertexSetNew.push_back(ProID);
			int rnew=rank[ProID],r;
			for(int i=0;i<ProBeginVertexSet.size();i++){
				r=rank[ProBeginVertexSet[i]];
				if(LCAQuery(rnew,r)!=rnew){
					ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
				}
			}
			ProBeginVertexSet=ProBeginVertexSetNew;
		}

	}

	int ProBeginVertexID;
	for(int i=0;i<ProBeginVertexSet.size();i++){
		ProBeginVertexID=ProBeginVertexSet[i];
		vector<int> linee; //linee.clear();
		linee.reserve(heightMax);
		int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
		while(Tree[rank[pachidd]].height>1){
			linee.insert(linee.begin(),pachidd);
			pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
		}
		linee.insert(linee.begin(),pachidd);

		eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum);
	}
	//return checknum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel){
	int childID=Tree[children].uniqueVertex;
	int childH=Tree[children].height-1;
	for(int i=0;i<Tree[children].dis.size();i++){
		if(Tree[children].cnt[i]==0){
			changelabel+=1;
			//firstly, check which dis can be infected
			int disBF=Tree[children].dis[i];
			int PID;
			//chidlID
			for(int k=0;k<VidtoTNid[childID].size();k++){
				PID=VidtoTNid[childID][k];
				if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
					Tree[PID].cnt[i]-=1;
				}
			}

			//line[i]
			for(int k=0;k<VidtoTNid[line[i]].size();k++){
				PID=VidtoTNid[line[i]][k];
				//if(PID>children){
//				if(Tree[PID].height>Tree[children].height){
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){
					if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
						Tree[PID].cnt[childH]-=1;
					}
				}
			}

			//secondly, calculate the actual distance
			int dis=INF; int count=0;
			int Dvb; int b,bH; int DDvb=INF;
			for(int j=0;j<Tree[children].vert.size();j++){
				Dvb=Tree[children].vert[j].second.first;
				b=Tree[children].vert[j].first;
				bH=Tree[rank[b]].height-1;
				if(bH<i){
					if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
						dis=Dvb+Tree[rank[line[i]]].dis[bH];
						count=1;
					}else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
						count+=1;
					}
				}else if(bH==i){
					DDvb=Dvb;
					if(Dvb<dis){
						dis=Dvb;
						count=1;
					}else if(Dvb==dis){
						count+=1;
					}
				}else{
					if(Dvb+Tree[rank[b]].dis[i]<dis){
						dis=Dvb+Tree[rank[b]].dis[i];
						count=1;
					}else if(Dvb+Tree[rank[b]].dis[i]==dis){
						count+=1;
					}
				}
			}
			if(DDvb==dis) Tree[children].FN[i]=true;
			Tree[children].dis[i]=dis;
			Tree[children].cnt[i]=count;
		}
	}

	line.push_back(childID);
	for(int i=0;i<Tree[children].ch.size();i++){
		eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel);
	}
	line.pop_back();
}
