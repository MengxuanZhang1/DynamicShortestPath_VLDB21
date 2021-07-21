/*
 * CHupd.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "head.h"

void Graph::NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
	sm->wait();
	int ID1, w1;
	int ID2, w2;
	for(int k=p.first;k<p.second;k++){
		ID1=Neighvec[k].first;
		for(int h=0;h<Neighvec.size();h++){
			ID2=Neighvec[h].first;
			insertEMTOrderGenerate(ID1, ID2, 1);
		}
	}
	sm->notify();
}

void Graph::insertEMTOrderGenerate(int u,int v,int w){
	if(E[u].find(v)==E[u].end()){
		E[u].insert(make_pair(v,make_pair(w,1)));
		DD[u]++;
		DD2[u]++;
	}
}

void Graph::NoRegularGraphOrder(){
	cout<<"begin order the degree!"<<endl;
	benchmark::heap<2, int, int> Q(nodenum);
	for(int i=0;i<nodenum;i++){
		Q.update(i,Neighbor[i].size());
		//cout<<"i: "<<i<<" Neighbor[i] size: "<<Neighbor[i].size()<<endl;
		if(Neighbor.size()==0) cout<<"Isolated node "<<i<<endl;
	}
	int topID, topdegree;
	int cnt=0;
	NodeOrder.assign(nodenum, 0);
	vNodeOrder.assign(nodenum, 0);
	while(!Q.empty()){
		Q.extract_min(topID, topdegree);
		NodeOrder[topID]=cnt;
		vNodeOrder[cnt]=topID;
		cnt+=1;
	}
	cout<<"Node finish ordering!"<<endl;
}

void Graph::ReadGraph(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}

	int nn,en;
	IF>>nodenum>>edgenum;

	eNum=0;
	set<pair<int,int>> eSet;

	vector<pair<int,int>> vecp; vecp.clear();
	Neighbor.assign(nodenum, vecp);

	set<int> setp; setp.clear();
	AdjacentNodes.assign(nodenum, setp);

	//to avoid the redundant information
	set<pair<int,int>> EdgeRedun;

	int ID1, ID2, weight;
	for(int i=0;i<edgenum;i++){
		IF>>ID1>>ID2>>weight;
		if(EdgeRedun.find(make_pair(ID1,ID2))==EdgeRedun.end()){
			Neighbor[ID1].push_back(make_pair(ID2, weight));
			AdjacentNodes[ID1].insert(ID2);
		}
		EdgeRedun.insert(make_pair(ID1,ID2));
	}
}

void Graph::ReadGraphForCHPinc(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}

	int nn,en;
	IF>>nodenum>>edgenum;

	eNum=0;
	set<pair<int,int>> eSet;

	vector<pair<int,int>> vecp; vecp.clear();
	Neighbor.assign(nodenum, vecp);

	//to avoid the redundant information
	set<pair<int,int>> EdgeRedun;

	int ID1, ID2, weight;
	for(int i=0;i<edgenum;i++){
		IF>>ID1>>ID2>>weight;
		if(EdgeRedun.find(make_pair(ID1,ID2))==EdgeRedun.end()){
			Neighbor[ID1].push_back(make_pair(ID2, weight));
		}
		EdgeRedun.insert(make_pair(ID1,ID2));

		if(ID1<ID2){
			if(eSet.find(make_pair(ID1,ID2))==eSet.end()){
				EdgeRe[make_pair(ID1,ID2)]=Edge.size();
				Edge.push_back(make_pair(ID1, ID2));
				eSet.insert(make_pair(ID1,ID2));
			}
		}
	}
	cout<<"Edge.size "<<Edge.size()<<" "<<nodenum<<" "<<edgenum<<endl;
	vSm.reserve(Edge.size());
	for(int i = 0; i < Edge.size(); i++)
		vSm[i] = new Semaphore(1);
}

void Graph::EffiCheckDij(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		//cout<<ID1<<" "<<ID2<<endl;
		ODpair.push_back(make_pair(ID1, ID2));
	}


	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	//for(int i=0;i<ODpair.size();i++){
	for(int i=0;i<200;i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=BiDij(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"BiDij query time "<<runT/ODpair.size()<<endl;
}

void Graph::EffiCheckCHPorder(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair.size();i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=QueryCHPorder(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CH query time "<<runT/ODpair.size()<<endl;
}

void Graph::EffiCheckCH(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair.size();i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=QueryCH(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CH query time "<<runT/ODpair.size()<<endl;
}

void Graph::EffiCheckH2H(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair.size();i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=QueryH2H(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"H2H query time "<<runT/ODpair.size()<<endl;
}

void Graph::EffiCheckPLL(string filename){
	ifstream IF(filename);
	if(!IF){
		cout<<"Cannot open Map "<<filename<<endl;
	}
	int num, ID1, ID2;
	vector<pair<int,int>> ODpair;
	IF>>num;
	for(int k=0;k<num;k++){
		IF>>ID1>>ID2;
		ODpair.push_back(make_pair(ID1, ID2));
	}

	int s, t, d;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<ODpair.size();i++){
		s=ODpair[i].first;
		t=ODpair[i].second;
		d=ShortestDisQuery(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL query time "<<runT/ODpair.size()<<endl;
}

void Graph::CorCheckDij(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=BiDij(s,t);
		if(d1!=d2)
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
	}
}

void Graph::CorCheckCHPorder(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=QueryCHPorder(s,t);
		if(d1!=d2)
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
	}
}

void Graph::CorCheckCH(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=QueryCH(s,t);
		if(d1!=d2)
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
	}
}

void Graph::CorCheckH2H(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=QueryH2H(s,t);
		if(d1!=d2)
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
	}
}

void Graph::CorCheckPLL(){
	srand (time(NULL));
	int s, t, d1, d2, d3;
	for(int i=0;i<100;i++){
		s=rand()%nodenum;
		t=rand()%nodenum;
		d1=Dij(s,t);
		d2=ShortestDisQuery(s,t);
		if(d1!=d2)
			cout<<"InCorrect!"<<s<<" "<<t<<" "<<d1<<" "<<d2<<endl;
	}
}

int Graph::Dij(int ID1, int ID2){
	if(ID1==ID2) return 0;
	//if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	benchmark::heap<2, int, int> pqueue(nodenum);
	pqueue.update(ID1,0);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);

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
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
				}
			}
		}
	}
	return d;
}

int Graph::BiDij(int ID1, int ID2){
	if(ID1==ID2) return 0;
	//if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;//to avoid the incorrectness caused by the isolated vertex
	benchmark::heap<2, int, int> queueF(nodenum), queueB(nodenum);
	queueF.update(ID1,0);
	queueB.update(ID2,0);

	vector<bool> closedF(nodenum, false), closedB(nodenum, false);
	vector<int> distanceF(nodenum, INF), distanceB(nodenum, INF);

	distanceF[ID1]=0;
	distanceB[ID2]=0;
	int topNodeIDF, topNodeDisF, topNodeIDB, topNodeDisB;
	int NNodeIDF,NWeighF, NNodeIDB, NWeighB;

	int d=INF;//initialize d to infinite for the unreachable case

	while(!queueF.empty() || !queueB.empty()){
		if(queueF.top()+queueB.top()>=d){
			return d;
		}
		//forward
		queueF.extract_min(topNodeIDF, topNodeDisF);
		closedF[topNodeIDF]=true;
		for(auto itF=Neighbor[topNodeIDF].begin();itF!=Neighbor[topNodeIDF].end();itF++){
			NNodeIDF=(*itF).first;
			NWeighF=(*itF).second+topNodeDisF;
			if(closedB[NNodeIDF] && NWeighF+distanceB[NNodeIDF]<d){
				d=NWeighF+distanceB[NNodeIDF];
			}
			if(!closedF[NNodeIDF]){
				if(distanceF[NNodeIDF]>NWeighF){
					distanceF[NNodeIDF]=NWeighF;
					queueF.update(NNodeIDF, NWeighF);
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
			}
			if(!closedB[NNodeIDB]){
				if(distanceB[NNodeIDB]>NWeighB){
					distanceB[NNodeIDB]=NWeighB;
					queueB.update(NNodeIDB, NWeighB);
				}
			}
		}
	}
	return d;
}

//A* algorithm is only implemented in Beijing Road Network
//Since only its coordinate is available
int Graph::Astar(int ID1, int ID2){
	if(ID1==ID2) return 0;
	benchmark::heap<2, int, int> pqueue(nodenum);
	int heurisDis=EuclideanDis(ID1,ID2);
	pqueue.update(ID1,heurisDis);

	vector<bool> closed(nodenum, false);
	vector<int> distance(nodenum, INF);//actual distance from source

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
				}
			}
		}
	}
	return d;
}

int Graph::EuclideanDis(int s, int t){
	int lat=(int)(abs(GraphLocation[s].first-GraphLocation[t].first)*111319);
	int lon=(int)(abs(GraphLocation[s].second-GraphLocation[t].second)*83907);
	int min,max;
	min=(lat>lon)?lon:lat;
	max=(lat>lon)?lat:lon;
	int approx=max*1007+min*441;
	if(max<(min<<4))
		approx-=max*40;
	return (approx+512)>>10;
}

int Graph::writeShortCutorder(string filename){
	//only write the one with higher order
	//write together with the higher adjacent vertex
	ofstream ofile(filename);
	cout<<"Writing shortcuts!"<<endl;
	for(int i = 0; i < nodenum; i++)
	{
		ofile << i << "\t" << NodeOrder[i] << "\t";

		vector<pair<int, int>> shortcut;
		int ID, Wei;
		for(int k=0;k<vvpShortCut[i].size();k++){
			ID=vvpShortCut[i][k].first;
			Wei=vvpShortCut[i][k].second;
			if(NodeOrder[ID]>NodeOrder[i])
				shortcut.push_back(make_pair(ID, Wei));
		}

		for(int j=0;j<Neighbor[i].size();j++){
			ID=Neighbor[i][j].first;
			Wei=Neighbor[i][j].second;
			if(NodeOrder[ID]>NodeOrder[i])
				shortcut.push_back(make_pair(ID, Wei));
		}

		if(shortcut.size() != 0)
		{
			ofile << 1;
			ofile << "\t" << shortcut.size();
			vector<pair<int,int>>::iterator iR;
			iR=shortcut.begin();
			while(iR!=shortcut.end()){
				ofile<<"\t"<<(*iR).first<<"\t"<<(*iR).second;
				iR++;
			}
		}
		else
			ofile << 0;

		ofile << endl;
	}
	return 0;
}

int Graph::ReadShortCut(string filename){
	ifstream inSCH(filename);
	int nodeID, nodeOrder, c, nsc, nrsc, n, d, m;
	cout<<"Reading shortcuts!"<<endl;

	for(int i = 0; i < nodenum; i++)
	{
		inSCH >> nodeID >> nodeOrder >> c;
		NodeOrder[nodeID] = nodeOrder;
		if(!c)
			continue;

		inSCH >> nrsc;
		for(int j = 0; j < nrsc; j++)
		{
			inSCH >> n >> d;
			AdjaShort[nodeID].push_back(make_pair(n,d));
			AdjaShortR[n].push_back(make_pair(nodeID,d));
		}
	}
	inSCH.close();

	cout<<"shortcut finish reading!"<<endl;
	return 0;
}

void Graph::writeCHPIncrease(string indexfile)
{
	ofstream ofSupportNodes(indexfile+"SupportNodes");
	ofstream ofPathInfor(indexfile+"PathInfor");
	ofstream ofEdgeOnPath(indexfile+"EdgeOnPath");

	for(int i = 0; i < nodenum; i++)
	{
		ofSupportNodes << i << "\t" << SupportNodes[i].size() << endl;
		for(auto& im : SupportNodes[i])
		{
			ofSupportNodes << im.first << "\t" << im.second.size();
			for(auto& iv : im.second)
				ofSupportNodes << "\t" << iv;
			ofSupportNodes << endl;
		}
	}
	ofSupportNodes.close();

	for(int i = 0; i < nodenum; i++)
	{
		ofPathInfor << i << "\t" << PathInfor[i].size() << endl;
		for(auto& im : PathInfor[i])
			ofPathInfor << im.first.first << "\t" << im.first.second << "\t" << im.second << endl;
	}
	ofPathInfor.close();

	for(int i = 0; i < (int)EdgeOnPath.size(); i++)
	{
		ofEdgeOnPath << i << "\t" << EdgeOnPath[i].size() << endl;
		for(auto& iTri : EdgeOnPath[i])
			ofEdgeOnPath << iTri.u << "\t" << iTri.v << "\t" << iTri.w << endl;
	}
	ofEdgeOnPath.close();
}

void Graph::readCHPIncrease(string indexfile)
{
	set<pair<int,int>> setpair;
	InvalidWP.assign(nodenum, setpair);
	ifstream ifSupportNodes(indexfile + "SupportNodes");
	ifstream ifPathInfor(indexfile + "PathInfor");
	ifstream ifEdgeOnPath(indexfile + "EdgeOnPath");

	int nID, ssize;
	for(int i = 0; i < nodenum; i++)
	{
		ifSupportNodes >> nID  >> ssize;
		int n, nsize;
		map<int, vector<int> > mv;
		for(int j = 0; j < ssize; j++)
		{
			ifSupportNodes >> n >> nsize;
			vector<int> vSupport(nsize, -1);
			for(int k = 0; k < nsize; k++)
				ifSupportNodes >> vSupport[k];
			mv[n] = vSupport;
		}

		SupportNodes[i] = mv;
	}
	ifSupportNodes.close();

	for(int i = 0; i < nodenum; i++)
	{
		ifPathInfor >> nID >> ssize;
		unordered_map<pair<int,int>, int, hash_pair> pathinfvec;
		int a, b, c;
		for(int j = 0; j < ssize; j++)
		{
			ifPathInfor >> a >> b >> c;
			pathinfvec.insert(make_pair(make_pair(a, b), c));
		}
		PathInfor[i] = pathinfvec;
	}
	ifPathInfor.close();

	int eID;
	for(int i = 0; i < (int)EdgeOnPath.size(); i++)
	{
		ifEdgeOnPath >> eID >> ssize;
		tri t;
		vector<tri> vt(ssize, t);
		for(int j = 0; j < ssize; j++)
		{
			tri t2;
			ifEdgeOnPath >> t2.u >> t2.v >> t2.w;
			vt[j] = t2;
		}
		EdgeOnPath[i] = vt;
	}
	ifEdgeOnPath.close();
}

int Graph::writeShortCutCH(string filename){

	ofstream ofile(filename);
	cout<<"Writing shortcuts!"<<endl;
	for(int i = 0; i < nodenum; i++)
	{
		ofile << i << "\t" << NodeOrder[i] << "\t";

		if(NeighborCon[i].size() != 0)
		{
			ofile << 1;

			ofile << "\t" << NeighborCon[i].size();
			vector<pair<int,pair<int,int>>>::iterator iR;
			iR=NeighborCon[i].begin();
			while(iR!=NeighborCon[i].end()){
				ofile<<"\t"<<(*iR).first<<"\t"<<(*iR).second.first<<"\t"<<(*iR).second.second;
				iR++;
			}
		}
		else
			ofile << 0;

		ofile << endl;
	}
	return 0;
}

int Graph::ReadShortCutCH(string filename){
	vector<pair<int,pair<int,int>>> vect;
	NeighborCon.assign(nodenum,vect);
	ifstream inSCH(filename);
	int nodeID, nodeOrder, c, nsc, nrsc, n, d, m;
	//cout<<"Reading shortcuts!"<<endl;

	for(int i = 0; i < nodenum; i++)
	{
		inSCH >> nodeID >> nodeOrder >> c;
		NodeOrder[nodeID] = nodeOrder;
		if(!c)
			continue;

		inSCH >> nrsc;
		for(int j = 0; j < nrsc; j++)
		{
			inSCH >> n >> d >> m;
			NeighborCon[nodeID].push_back(make_pair(n,make_pair(d,m)));
		}
	}
	inSCH.close();

	cout<<"shortcut finish reading!"<<endl;
	return 0;
}

void Graph::writePLL(string filename, string filenameP){
	ofstream OF(filename);
	OF<<nodenum<<endl;
	for(int nodeID=0;nodeID<nodenum;nodeID++){
		OF<<nodeID<<" "<<NodeOrder[nodeID]<<" "<<Label[nodeID].size();
		for(unordered_map<int,int>::iterator it=Label[nodeID].begin();it!=Label[nodeID].end();it++){
			OF<<" "<<(*it).first<<" "<<(*it).second;
		}
		OF<<endl;
	}

	ofstream OF1(filenameP);
	for(int nodeID=0;nodeID<nodenum;nodeID++){
		OF1<<nodeID<<" "<<PruningPointNew[nodeID].size();
		for(auto itp=PruningPointNew[nodeID].begin();itp!=PruningPointNew[nodeID].end();itp++){
			OF1<<" "<<(*itp).first<<" "<<(*itp).second.size();
			for(int k=0;k<(*itp).second.size();k++){
				OF1<<" "<<(*itp).second[k];
			}
		}
		OF1<<endl;
	}

	//cout<<"PLL index finish writing!"<<endl;
}

void Graph::readPLL(string filename, string filenameP){
	ifstream IF(filename);
	int num;
	IF>>num;
	unordered_map<int,int> M0;
	Label.assign(num,M0);
	NodeOrder.assign(num,0);
	vNodeOrder.assign(num,0);

	int nodeID, order, pnum, hubid, dis;
	for(int i=0;i<num;i++){
		IF>>nodeID>>order>>pnum;
		NodeOrder[nodeID]=order;
		vNodeOrder[order]=nodeID;
		for(int j=0;j<pnum;j++){
			IF>>hubid>>dis;
			Label[nodeID].insert(make_pair(hubid,dis));
		}
	}

	unordered_map<int,vector<int>> unorderm;
	PruningPointNew.assign(num,unorderm);
	ifstream IF1(filenameP);

	int pairnum, pairsize;
	int c,u;
	for(int i=0;i<num;i++){
		IF1>>nodeID>>pairnum;
		for(int p=0;p<pairnum;p++){
			IF1>>c>>pairsize;
			vector<int> vec;
			for(int k=0;k<pairsize;k++){
				IF1>>u;
				vec.push_back(u);
			}
			PruningPointNew[nodeID][c]=vec;
		}
	}

	//cout<<"PLL index finish reading!"<<endl;
}

void Graph::ReadIns(string filename,vector<pair<pair<int,int>,int>>& TestData){
	TestData.clear();

	int num, ID1, ID2, neww;
	ifstream IF(filename);
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>ID1>>ID2>>neww;
		TestData.push_back(make_pair(make_pair(ID1, ID2), neww));
	}
	IF.close();
}

/*Index Size Computation*/
long long Graph::H2Hindexsize(){
	long long m=0;

	for(int i=0;i<Tree.size();i++){
		m+=Tree[i].dis.size()*sizeof(int);
	}

	cout<<"H2H index size "<<(double)m/1024/1024<<endl;
	return m;
};

long long Graph::CHindexsize(){
	long long m=0;

	for(int i=0;i<NeighborCon.size();i++){
		m+=NeighborCon[i].size()*3*sizeof(int);
	}

	cout<<"CH index size "<<(double)m/1024/1024<<endl;
	return m;
}

//SCconNodesMT size for both CH and H2H increase
long long Graph::SCconNodesize(){
	long long m=0;

	for(int i=0;i< SCconNodesMT.size();i++){
		for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
			m+=sizeof(int)+(*it).second.size()*sizeof(int);
		}
	}

	cout<<"CH Support Nodes size "<<(double)m/1024/1024<<endl;
	return m;
}

long long Graph::CHPindexsize(){
	long long m=0;

	for(int i=0;i<AdjaShort.size();i++){
		m+=AdjaShort[i].size()*2*sizeof(int);
	}

	cout<<"CHP index size "<<(double)m/1024/1024<<endl;
	return m;
}

//witness path size for CHP increase
long long Graph::CHPWitPathsize(){
	long long m1=0,m2=0;

	for(int i=0;i< SupportNodes.size();i++){
		for(auto it=SupportNodes[i].begin(); it!=SupportNodes[i].end(); it++){
			m1+=sizeof(int)+(*it).second.size()*sizeof(int);
		}
	}
	cout<<"CHP Support Nodes size "<<(double)m1/1024/1024<<endl;

	for(int i=0;i<PathInfor.size();i++){
		for(auto it=PathInfor[i].begin(); it!=PathInfor[i].end(); it++){
			m2+=sizeof(int)*3;
		}
	}

	for(int i=0;i<EdgeOnPath.size();i++){
		m2+=EdgeOnPath[i].size()*3*sizeof(int);
	}

	cout<<"CHP witness path size "<<(double)m2/1024/1024<<endl;
	cout<<"Support+Witness size "<<(double)(m1+m2)/1024/1024<<endl;
	return 0;
}

long long Graph::PLLindexsize(){
	long long m=0;
	for(int nodeID=0;nodeID<nodenum;nodeID++){
		m+=Label[nodeID].size()*2*sizeof(int);
	}
	cout<<"PLL index size "<<(double)m/1024/1024<<endl;
	return m;
}

//pruning point size for PLL-P increase
long long Graph::PrunePointsize(){
	long long m=0;
	for(int nodeID=0;nodeID<nodenum;nodeID++){
		for(auto it=PruningPointNew[nodeID].begin(); it!=PruningPointNew[nodeID].end(); it++){
			m+=sizeof(int)+(*it).second.size()*sizeof(int);
		}
	}
	cout<<"PLL Pruning Point size "<<(double)m/1024/1024<<endl;
	return m;
}


int	Graph::QueryCHPorder(int ID1, int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

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
			//cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
					//cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
				}
			}

			for(auto out=AdjaShort[topNodeIDForward].begin();out!=AdjaShort[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						//if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
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
					//cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
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
					}
				}
			}
		}
	}
	return d;
}

int	Graph::QueryCH(int ID1, int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	benchmark::heap<2,int,int> fHeapForward(nodenum);
	benchmark::heap<2, int, int> fHeapBackward(nodenum);

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
			//cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

			if(vDistanceForward[topNodeIDForward] > d)
				bF = true;

			vVisitedF[topNodeIDForward] = true;

			if(vVisitedB[topNodeIDForward]){
				int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
				if(distTmp<d){
					d=distTmp;
					//cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
				}
			}

			for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
				neighborNodeID = (*out).first;
				neighborLength = (*out).second.first;

				int df = vDistanceForward[topNodeIDForward] + neighborLength;
				if(!vVisitedF[neighborNodeID]){
					if(vDistanceForward[neighborNodeID] > df){
						//if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
						vDistanceForward[neighborNodeID] = df;
						fHeapForward.update(neighborNodeID, df);
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
					//cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
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
					}
				}
			}
		}
	}
	return d;
}

int Graph::QueryH2H(int ID1,int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int r1=rank[ID1], r2=rank[ID2];
	int LCA=LCAQuery(r1,r2);

	if(LCA==r1)
		return Tree[r2].dis[Tree[r1].pos.back()];
	else if(LCA==r2)
		return Tree[r1].dis[Tree[r2].pos.back()];
	else{
		int tmp=INF;
		for(int i=0;i<Tree[LCA].pos.size();i++){
			if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
				tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
		}
		return tmp;
	}
}

int Graph::ShortestDisQuery(int ID1,int ID2){
	if(ID1==ID2) return 0;
	if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
	int d=INF;
	unordered_map<int,int>::iterator it;
	int hub, dis1, dis2;
	for(it=Label[ID1].begin();it!=Label[ID1].end();it++){
		hub=(*it).first;
		dis1=(*it).second;
		if(Label[ID2].find(hub)!=Label[ID2].end()){
			dis2=Label[ID2][hub];
			if(dis1+dis2<d){
				d=dis1+dis2;
				//cout<<"details "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
			}
		}
	}
	return d;
}
