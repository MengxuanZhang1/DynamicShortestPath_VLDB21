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
	string orderfile=DesFile+"Order";
	g.ReadGraph(destiGraph);

	//index construction
	t1=std::chrono::high_resolution_clock::now();
	g.CHconsorderMT(orderfile);
//    g.CHconsMTOrderGenerate(orderfile);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CH construction time "<<runT<<endl;
    exit(0);
	g.CorCheckCH();

	//index size & efficiency
	string ODfile=DesFile+"OD";
	g.CHindexsize();
	g.SCconNodesize();
	g.EffiCheckCH(ODfile);

	//index update data
	string updateFile=DesFile+"Update";
	vector<pair<pair<int,int>,int>> testdata;
	g.ReadIns(updateFile, testdata);
	vector<pair<pair<int,int>,pair<int,int>>> testdataInc, testdataDec;
	for(int k=0;k<testdata.size();k++){
		testdataInc.push_back(make_pair(testdata[k].first,make_pair(testdata[k].second,testdata[k].second*10)));
		testdataDec.push_back(make_pair(testdata[k].first,make_pair(testdata[k].second,testdata[k].second*0.5)));
	}

	//CHW decrease
	Graph g1=g;
	t1=std::chrono::high_resolution_clock::now();
	g1.CHdecBat(testdataDec);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"DCH decrease "<<runT<<" "<<runT/testdata.size()<<endl;

	//CHW increase
	Graph g2=g;
	t1=std::chrono::high_resolution_clock::now();
	g2.CHincBatMT(testdataInc);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"DCH increase "<<runT<<" "<<runT/testdata.size()<<endl;

	//UE decrease
	Graph g3=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g3.CHdecStr(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*0.5);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"UE decrease "<<runT<<" "<<runT/testdata.size()<<endl;

	//UE increase
	Graph g4=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g4.CHincStrMT(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*10);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"UE increase "<<runT<<" "<<runT/testdata.size()<<endl;

	return 0;
}



