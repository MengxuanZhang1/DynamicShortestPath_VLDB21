/*
 * Labelcon.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "head.h"

int Graph::StartPLL(string filename, string filenameP, string orderfile){
	fstream file;
	file.open(filename);
	if(!file){
		std::chrono::high_resolution_clock::time_point t1, t2;
		std::chrono::duration<double> time_span;
		double runT;

		t1=std::chrono::high_resolution_clock::now();
		PLLcon(orderfile);
		t2=std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
		runT= time_span.count();
		cout<<"PLL construction time "<<runT<<endl;

		writePLL(filename, filenameP);
		readPLL(filename, filenameP);
	}
	else
	{
		readPLL(filename, filenameP);
	}

	return 0;
}

void Graph::PLLcon(string orderfile){
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
	Label.assign(nodenum,m);
	unordered_map<int,vector<int>> unorderm;
	PruningPointNew.assign(nodenum,unorderm);

	int ID;
	int cnt=0;
	for(int i=nodenum-1;i>=0;i--){
		ID=vNodeOrder[i];
		vector<pair<int,int>> vp;
		DijksPrune1(ID,vp);
		cnt+=1;
		for(int j=0;j<vp.size();j++){
			Label[vp[j].first].insert(make_pair(ID, vp[j].second));
			//cout<<vp[j].first<<" "<<vp[j].second<<endl;
		}
	}
}

void Graph::DijksPrune1(int nodeID, vector<pair<int,int>>& vp){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(nodeID,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

	distance[nodeID]=0;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		closed[topNodeID]=true;

		int TempDis; vector<int> SupNode;
		ShortestDisQuery1(nodeID, topNodeID,SupNode,TempDis);
		if(TempDis<=topNodeDis){

			if(NodeOrder[topNodeID]<NodeOrder[nodeID]){
				for(int k=0;k<SupNode.size();k++){
					int supn=SupNode[k];

					PruningPointNew[topNodeID][supn].push_back(nodeID);
					PruningPointNew[nodeID][supn].push_back(topNodeID);
				}
			}
			continue;
		}


		//Label[topNodeID].insert(nodeID, topNodeDis);
		vp.push_back(make_pair(topNodeID,topNodeDis));
		for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				}
			}
		}
	}
}

void Graph::PLLdec(int a,int b, int oldW, int newW){
	for(int i=0;i<Neighbor[a].size();i++){
		if(Neighbor[a][i].first==b){
			Neighbor[a][i].second=newW;
			break;
		}
	}
	for(int j=0;j<Neighbor[b].size();j++){
		if(Neighbor[b][j].first==a){
			Neighbor[b][j].second=newW;
			break;
		}
	}

	benchmark::heap<2, int, int> Q0(nodenum);
	unordered_map<int,int>::iterator it1,it2;
	int hubID;
	for(it1=Label[a].begin();it1!=Label[a].end();it1++){
		hubID=(*it1).first;
		Q0.update(hubID, nodenum-NodeOrder[hubID]);
		//Q0.update(hubID, NodeOrder[hubID]);
	}
	for(it2!=Label[b].begin();it2!=Label[b].end();it2++){
		hubID=(*it2).first;
		Q0.update(hubID, nodenum-NodeOrder[hubID]);
		//Q0.update(hubID, NodeOrder[hubID]);
	}

	int NodeID, NodeR;
	while(!Q0.empty()){
		Q0.extract_min(NodeID, NodeR);
		if(Label[a].find(NodeID)!=Label[a].end()){
			DijkstraPrune1(NodeID, b, Label[a][NodeID]+newW);
		}
		if(Label[b].find(NodeID)!=Label[b].end()){
			DijkstraPrune1(NodeID, a, Label[b][NodeID]+newW);
		}
	}
}

void Graph::PLLinc(int u,int v, int oldW, int newW){
	vector<int> Disu, Disv;
	vector<int> DisuNew, DisvNew;
	DijkstraList(u, Disu);
	DijkstraList(v, Disv);

	for(int i=0;i<Neighbor[u].size();i++){
		if(Neighbor[u][i].first==v){
			Neighbor[u][i].second=newW;
			break;
		}
	}
	for(int j=0;j<Neighbor[v].size();j++){
		if(Neighbor[v][j].first==u){
			Neighbor[v][j].second=newW;
			break;
		}
	}

	DijkstraList(u, DisuNew);
	DijkstraList(v, DisvNew);

	//identify the possible & real affected labels
	vector<int> PAu, PAv, RAu, RAv;
	PLLaffect(u, oldW, PAu, RAu, Disu, Disv, DisvNew);
	PLLaffect(v, oldW, PAv, RAv, Disv, Disu, DisuNew);

	//identify the invalid and missing labels
	vector<int> ILu,ILv,MLu,MLv;
	PLLinvalid(u, ILu, MLu, RAu, RAv, PAv);
	PLLinvalid(v, ILv, MLv, RAv, RAu, PAu);

	//remove the invalid label
	PLLremove(u, RAu, ILu);
	PLLremove(v, RAv, ILv);

	//connected or not after edge deletion
	int disuv=Dij(u,v);
	if(disuv!=INF){//still connected
		//collect IL&ML in decreasing vertex order
		//collect RA in a set
		set<int> Affect;
		for(int i=0;i<RAu.size();i++)
			Affect.insert(RAu[i]);
		for(int i=0;i<RAv.size();i++)
			Affect.insert(RAv[i]);

		map<int,int> Rmap;
		for(int i=0;i<ILu.size();i++)
			Rmap.insert(make_pair(NodeOrder[ILu[i]],ILu[i]));
		for(int i=0;i<ILv.size();i++)
			Rmap.insert(make_pair(NodeOrder[ILv[i]],ILv[i]));
		for(int i=0;i<MLu.size();i++)
			Rmap.insert(make_pair(NodeOrder[MLu[i]],MLu[i]));
		for(int i=0;i<MLv.size();i++)
			Rmap.insert(make_pair(NodeOrder[MLv[i]],MLv[i]));

		int rOrder, rID;
		while(Rmap.size()>0){
			rOrder=(*Rmap.rbegin()).first;
			rID=(*Rmap.rbegin()).second;
			PLLPrune(rID, Affect);
			Rmap.erase(rOrder);
		}
	}
}

void Graph::PLLPrune(int r, set<int> Affect){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(r,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

	distance[r]=0;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		closed[topNodeID]=true;

		int TempDis=ShortestDisQuery(r, topNodeID);

		if(TempDis>topNodeDis){
			/*if(Affect.find(topNodeID)!=Affect.end()){
				if(Label[topNodeID][r]!=topNodeDis)
					cout<<"inconsistent///////////////////////////"<<endl;
			}*/
			Label[topNodeID][r]=topNodeDis;
			for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
				NNodeID=(*it).first;
				NWeigh=(*it).second+topNodeDis;
				if(!closed[NNodeID]){
					if(distance[NNodeID]>NWeigh){
						distance[NNodeID]=NWeigh;
						pqueue.update(NNodeID, NWeigh);
					}
				}
			}
		}

	}
}

void Graph::PLLremove(int u, vector<int> RAu, vector<int> ILu){
	int s,r;
	set<int> ILuMap;
	for(int k=0;k<ILu.size();k++)
		ILuMap.insert(ILu[k]);

	for(int i=0;i<RAu.size();i++){
		s=RAu[i];
		vector<int> R;
		unordered_map<int,int>::iterator it;
		for(it=Label[s].begin();it!=Label[s].end();it++){
			r=(*it).first;
			if(ILuMap.find(r)!=ILuMap.end()){
				R.push_back(r);
			}
		}
		for(int j=0;j<R.size();j++){
			Label[s].erase(R[j]);
		}
	}
}

void Graph::PLLinvalid(int u, vector<int>& ILu, vector<int>& MLu, vector<int> RAu, vector<int> RAv, vector<int> PAv){
	int t,r,w,ID;
	for(int i=0;i<RAu.size();i++){
		t=RAu[i];
		for(int j=0;j<RAv.size();j++){
			r=RAv[j];
			if(Label[t].find(r)!=Label[t].end()){
				ILu.push_back(r);

				//find ri
				int order=0;
				int ri=-1;
				unordered_map<int,int>::iterator it;
				for(it=Label[t].begin();it!=Label[t].end();it++){
					ID=(*it).first;
					if(NodeOrder[ID]<NodeOrder[r]){
						if(NodeOrder[ID]>order){
							order=NodeOrder[ID];
							ri=ID;
						}
					}
				}
				if(ri!=-1){
					set<int> PAvmap;
					for(int k=0;k<PAv.size();k++)
						PAvmap.insert(PAv[k]);

					for(int Order=NodeOrder[ri]+1;Order<NodeOrder[r];Order++){
						w=vNodeOrder[Order];
						if(PAvmap.find(w)!=PAvmap.end()){
							int disWT=ShortestDisQuery(w,t);
							if(Label[w].find(r)!=Label[w].end() && Label[w][r]+Label[t][r]==disWT){
								MLu.push_back(w);
							}
						}
					}
				}

			}
		}
	}

}

void Graph::PLLaffect(int u, int oldw, vector<int>& PAu, vector<int>& RAu, vector<int> disu, vector<int> disv, vector<int> disvNew){
	vector<bool> Flag(nodenum, false);
	Flag[u]=true;
	queue<int> Q;
	Q.push(u);
	int top,r;
	while(!Q.empty()){
		top=Q.front();
		Q.pop();
		for(int nid=0; nid<Neighbor[top].size(); nid++){
			r=Neighbor[top][nid].first;
			//cout<<"r "<<r<<endl;
			if(!Flag[r]){
				if(disv[r]==disu[r]+oldw){
					Q.push(r);
					PAu.push_back(r);
					//cout<<"PA push"<<endl;
					if(disvNew[r]!=disu[r]+oldw){
						RAu.push_back(r);
					}
					Flag[r]=true;
				}
			}
		}
	}
}

void Graph::PSLdec(int a,int b, int oldW, int newW){
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

	//Chnum=0;
	//check the dis(a,b)
	int Dab=ShortestDisQuery(a,b);

	int LID,HID;
	if(NodeOrder[a]>NodeOrder[b]){
		LID=b; HID=a;
	}else{
		LID=a; HID=b;
	}

	if(Dab>newW){//the index change is triggered
		vector<vector<pair<int,int>>> Change;
		vector<pair<int,int>> vec;
		Change.assign(nodenum,vec);
		set<int> WaitPro;

		Label[LID][HID]=newW;
		//Chnum+=1;
		Change[LID].push_back(make_pair(HID,newW));
		WaitPro.insert(LID);

		//check the label of a,b
		int hubid, hubdis;
		unordered_map<int,int>::iterator it1=Label[LID].begin();
		for(;it1!=Label[LID].end();it1++){
			hubid=(*it1).first; hubdis=(*it1).second;
			if(NodeOrder[hubid]>NodeOrder[HID] && newW+hubdis<ShortestDisQuery(HID,hubid)){
				Label[HID][hubid]=newW+hubdis;
				//Chnum+=1;
				Change[HID].push_back(make_pair(hubid, newW+hubdis));
				WaitPro.insert(HID);
			}
		}
		unordered_map<int,int>::iterator it2=Label[HID].begin();
		for(;it2!=Label[HID].end();it2++){
			hubid=(*it2).first; hubdis=(*it2).second;
			if(newW+hubdis<ShortestDisQuery(LID, hubid)){
				Label[LID][hubid]=newW+hubdis;
				//Chnum+=1;
				Change[LID].push_back(make_pair(hubid, newW+hubdis));
				WaitPro.insert(LID);
			}
		}

		//check the label of their neighbors step by step
		while(WaitPro.size()>0){
			set<int> WaitProTem;
			vector<vector<pair<int,int>>> ChangeTem;
			vector<pair<int,int>> vec;
			ChangeTem.assign(nodenum,vec);

			for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
				int curID=*it;
				vector<pair<int,int>> curChange=Change[curID];
				int neiID, neiDis, hID, hDis;
				for(int i=0;i<Neighbor[curID].size();i++){
					neiID=Neighbor[curID][i].first;
					neiDis=Neighbor[curID][i].second;

					for(int j=0;j<curChange.size();j++){
						hID=curChange[j].first; hDis=curChange[j].second;
						if(NodeOrder[hID]>NodeOrder[neiID] && ShortestDisQuery(neiID, hID)>neiDis+hDis){
							Label[neiID][hID]=neiDis+hDis;
							//Chnum+=1;
							WaitProTem.insert(neiID);
							ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
						}
					}
				}
			}

			WaitPro=WaitProTem;
			Change=ChangeTem;
		}
	}
}

void Graph::PSLinc(int a,int b, int oldW, int newW){
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

	//Chnum=0;

	int LID,HID;
	if(NodeOrder[a]>NodeOrder[b]){
		LID=b; HID=a;
	}else{
		LID=a; HID=b;
	}

	int dis,disvally,dispeak;
	//activate or not
	dis=DisQueryLower1(LID,HID);

	//if(dispeak<=oldW) return;
	if(Label[LID].find(HID)!=Label[LID].end() && dis>oldW){//index update triggered
		vector<vector<pair<int,int>>> Change;
		vector<pair<int,int>> vec;
		Change.assign(nodenum,vec);
		set<int> WaitPro;
		vector<vector<int>> ChangeP;
		vector<int> vecint;
		ChangeP.assign(nodenum,vecint);
		set<int> WaitProP;

		WaitPro.insert(LID);
		Change[LID].push_back(make_pair(HID, oldW));
		disvally=DisQueryVally(LID,HID);
		Label[LID][HID]=disvally;//correct to the new value
		//Chnum+=1;

		//affected by the w(a.b)
		int hubID, hDis;
		int dis, cnt;
		for(unordered_map<int,int>::iterator it=Label[HID].begin();it!=Label[HID].end();it++){
			hubID=(*it).first; hDis=(*it).second;
			if(hubID!=HID && Label[LID].find(hubID)!=Label[LID].end() && oldW+hDis==Label[LID][hubID]){
				disvally=DisQueryVally(LID,hubID);
				if(Label[LID][hubID]<disvally){
					WaitPro.insert(LID);
					Change[LID].push_back(make_pair(hubID, oldW+hDis));
					Label[LID][hubID]=disvally;
					//Chnum+=1;
				}
			}
		}

		for(unordered_map<int,int>::iterator it=Label[LID].begin();it!=Label[LID].end();it++){
			hubID=(*it).first; hDis=(*it).second;
			if(Label[HID].find(hubID)!=Label[HID].end() && oldW+hDis==Label[HID][hubID]){
				disvally=DisQueryVally(HID,hubID);
				if(Label[HID][hubID]<disvally){
					WaitPro.insert(HID);
					Change[HID].push_back(make_pair(hubID, oldW+hDis));
					Label[HID][hubID]=disvally;
					//Chnum+=1;
				}
			}
		}

		while(WaitProP.size()>0 || WaitPro.size()>0){
			set<int> WaitProTem;
			vector<vector<pair<int,int>>> ChangeTem;
			vector<pair<int,int>> vec;
			ChangeTem.assign(nodenum,vec);
			set<int> WaitProPTem;
			vector<int> vecint;
			vector<vector<int>> ChangePTem;
			ChangePTem.assign(nodenum,vecint);

			//Change->Change & ChangeP
			for(set<int>::iterator it=WaitPro.begin();it!=WaitPro.end();it++){
				int curID=*it;
				vector<pair<int,int>> curChange=Change[curID];
				int neiID, neiDis, hID, hDis;

				for(int j=0;j<curChange.size();j++){
					hID=curChange[j].first; hDis=curChange[j].second;
					//Change->Change
					for(int k=0;k<Neighbor[curID].size();k++){
						neiID=Neighbor[curID][k].first; neiDis=Neighbor[curID][k].second;
						if(Label[neiID].find(hID)!=Label[neiID].end() && neiDis+hDis==Label[neiID][hID]){
							disvally=DisQueryVally(neiID,hID);
							if(Label[neiID][hID]<disvally){
								WaitProTem.insert(neiID);
								ChangeTem[neiID].push_back(make_pair(hID, neiDis+hDis));
								Label[neiID][hID]=disvally;
								//Chnum+=1;
							}
						}
					}

					//Change->ChangeP
					unordered_map<int,vector<int>>::iterator itt;
					if((itt=PruningPointNew[curID].find(hID))!=PruningPointNew[curID].end()){
						vector<int> vec=(*itt).second;
						for(int snum=0;snum<vec.size();snum++){
							int s=vec[snum];
					//if(PruningPoint.find(make_pair(curID,hID))!=PruningPoint.end()){
						//for(int snum=0;snum<PruningPoint[make_pair(curID,hID)].size();snum++){
							//int s=PruningPoint[make_pair(curID,hID)][snum];
							if(NodeOrder[s]<NodeOrder[curID] && NoSupportedPair.find(make_pair(s,curID))==NoSupportedPair.end()){
								disvally=DisQueryVally(s,curID);
								dispeak=DisQueryPeak(s,curID);
								if(dispeak>disvally){
									WaitProPTem.insert(s);
									ChangePTem[s].push_back(curID);
									Label[s][curID]=disvally;
									//Chnum+=1;
									NoSupportedPair.insert(make_pair(s,curID));
								}
							}else if(NodeOrder[s]>NodeOrder[curID] && NoSupportedPair.find(make_pair(curID,s))==NoSupportedPair.end()){
								disvally=DisQueryVally(curID,s);
								dispeak=DisQueryPeak(curID,s);
								if(dispeak>disvally){
									WaitProPTem.insert(curID);
									ChangePTem[curID].push_back(s);
									Label[curID][s]=disvally;
									//Chnum+=1;
									NoSupportedPair.insert(make_pair(curID,s));
								}
							}
						}
					}
				}
			}

			//ChangeP->CHangeP
			int v,u,neiid,neiw;
			for(set<int>::iterator itp=WaitProP.begin();itp!=WaitProP.end();itp++){
				v=*itp;
				for(int k=0;k<ChangeP[v].size();k++){
					u=ChangeP[v][k];
					for(int l=0;l<Neighbor[v].size();l++){
						neiid=Neighbor[v][l].first; neiw=Neighbor[v][l].second;
						if(NodeOrder[neiid]<NodeOrder[u]){
							disvally=DisQueryVally(neiid, u);
							dispeak=DisQueryPeak(neiid, u);
							if(disvally<dispeak){
								if(Label[neiid].find(u)==Label[neiid].end() || (Label[neiid].find(u)!=Label[neiid].end() && Label[neiid][u]!=disvally)){
									WaitProPTem.insert(neiid);
									ChangePTem[neiid].push_back(u);
									Label[neiid][u]=disvally;
									//Chnum+=1;
									NoSupportedPair.insert(make_pair(neiid,u));
								}
							}
						}
					}

				}
			}

			WaitPro=WaitProTem;
			Change=ChangeTem;
			WaitProP=WaitProPTem;
			ChangeP=ChangePTem;
		}

	}
}

int Graph::DijkstraList(int ID1, vector<int>& distance){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	distance.assign(nodenum, INF);
//	vector<int> prece(nodenum, 0);
	distance[ID1]=0;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		closed[topNodeID]=true;

		for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				//	prece[NNodeID]=topNodeID;
				}
			}
		}
	}
	return 0;
}

void Graph::DijkstraPrune1(int NodeID, int ID, int d){
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID,d);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

	distance[ID]=d;
	int topNodeID, topNodeDis;
	vector<pair<int,int>>::iterator it;
	int NNodeID,NWeigh;

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		closed[topNodeID]=true;

		int TempDis=PrefixalDisQuery(NodeID, topNodeID);
		if(TempDis<=topNodeDis)
			continue;

		Label[topNodeID][NodeID]=topNodeDis;
		for(it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				}
			}
		}
	}
}

int Graph::ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d){
	d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(Label[ID2].find(hub)!=Label[ID2].end()){
			dis2=Label[ID2][hub];
			if(dis1+dis2<d){
				d=dis1+dis2;
				SupNode.clear();
				SupNode.push_back(hub);
			}else if(dis1+dis2==d){
				SupNode.push_back(hub);
			}
		}
	}

	return d;
}

int Graph::DisQueryVally(int ID1, int ID2){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<Neighbor[ID1].size();i++){
		neiID=Neighbor[ID1][i].first;
		neiDis=Neighbor[ID1][i].second;
		if(NodeOrder[neiID]<=NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){
			if(neiDis+Label[neiID][ID2]<d){
				d=neiDis+Label[neiID][ID2];
			}
		}
	}
	return d;
}

int Graph::DisQueryLower1(int ID1, int ID2){
	int neiID,neiDis;
	int d=INF;
	for(int i=0;i<Neighbor[ID1].size();i++){
		neiID=Neighbor[ID1][i].first;
		neiDis=Neighbor[ID1][i].second;
		if(NodeOrder[neiID]<NodeOrder[ID2] && Label[neiID].find(ID2)!=Label[neiID].end()){
			if(neiDis+Label[neiID][ID2]<d){
				d=neiDis+Label[neiID][ID2];
			}
		}
	}
	return d;
}

int Graph::DisQueryPeak(int ID1, int ID2){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(hub!=ID2 && Label[ID2].find(hub)!=Label[ID2].end()){
			if(dis1+Label[ID2][hub]<d){
				d=dis1+Label[ID2][hub];
			}
		}
	}
	return d;
}

int Graph::PrefixalDisQuery(int ID1, int ID2){
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(NodeOrder[hub]<=NodeOrder[ID1] && Label[ID2].find(hub)!=Label[ID2].end()){
			dis2=Label[ID2][hub];
			if(dis1+dis2<d){
				d=dis1+dis2;
				//cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
			}
		}
	}
	return d;
}
