/*
 * CHupd.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "head.h"

int main(){
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

	string DesFile="./data/";
	vector<pair<pair<int,int>,int>> Testdata;

	Graph g;
	//read graph
	string destiGraph=DesFile+"Graphsw10k5";
	g.ReadGraph(destiGraph);
	g.CorCheckDij();
	string ODfile=DesFile+"OD";
	g.EffiCheckDij(ODfile);

	return 0;
}


