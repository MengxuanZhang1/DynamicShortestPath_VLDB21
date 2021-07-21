/*
 * CHupd.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "head.h"

int main(){
	string DesFile="./data/";
	string destiGraph=DesFile+"Graphsw10k5";
	string orderfile=DesFile+"Order";
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	Graph g;
	g.ReadGraph(destiGraph);
	srand (time(NULL));

	//time test for Bi-Dijkstra/A* shortest path query
	//firstly read the order file
	ifstream IF(orderfile);
	if(!IF)
		cout<<"Cannot open Map "<<orderfile<<endl;
	g.NodeOrder.assign(g.nodenum, -1);
	g.vNodeOrder.assign(g.nodenum, -1);
	int num, nodeID, nodeorder;
	IF>>num;
	for(int i=0;i<num;i++){
		IF>>nodeID>>nodeorder;
		g.NodeOrder[nodeID]=nodeorder;
		if(nodeorder!=-1){
			g.vNodeOrder[nodeorder]=nodeID;
		}
	}
	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<1000;i++){
		int s=rand()%g.nodenum;
		int t=rand()%g.nodenum;
		g.BiDijPath(s,t);
		//g.AstarPath(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Path retrieve time "<<runT<<" "<<runT/1000<<endl;

	//time test for PLL shortest path query
	string indexfile=DesFile+"/PLLPathIndex";
	g.PLLconsPath(orderfile,indexfile);
	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<1000;i++){
		int s=rand()%g.nodenum;
		int t=rand()%g.nodenum;
		g.PathRetriPLL(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Path retrieve time "<<runT<<" "<<runT/1000<<endl;

	//time test for CHP shortest path query
	g.StartCHPPath(orderfile);
	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<1000;i++){
		int s=rand()%g.nodenum;
		int t=rand()%g.nodenum;
		g.PathRetriCHP(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Path retrieve time "<<runT<<" "<<runT/1000<<endl;

	//time test for CHW shortest path query
	g.CHconsPath(orderfile);
	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<1000;i++){
		int s=rand()%g.nodenum;
		int t=rand()%g.nodenum;
		g.PathRetriCH(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Path retrieve time "<<runT<<" "<<runT/1000<<endl;

	//time test for H2H shortest path query
	g.H2HconPath(orderfile);
	t1=std::chrono::high_resolution_clock::now();
	for(int i=0;i<1000;i++){
		int s=rand()%g.nodenum;
		int t=rand()%g.nodenum;
		g.PathRetriH2H(s,t);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Path retrieve time "<<runT<<" "<<runT/1000<<endl;

	return 0;
}



