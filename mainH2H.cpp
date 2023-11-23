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
	g.H2HconOrderMT(orderfile);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"H2H construction time "<<runT<<endl;
	g.CorCheckH2H(10000);

	//index size & efficiency
	string ODfile=DesFile+"OD";
	g.H2Hindexsize();
	g.EffiCheckH2H(ODfile);

	//index update data
	string updateFile=DesFile+"Update";
	vector<pair<pair<int,int>,int>> testdata;
	g.ReadIns(updateFile, testdata);
	vector<pair<pair<int,int>,pair<int,int>>> testdataInc, testdataDec;
	for(int k=0;k<testdata.size();k++){
		testdataInc.push_back(make_pair(testdata[k].first,make_pair(testdata[k].second,testdata[k].second*1.5)));
		testdataDec.push_back(make_pair(testdata[k].first,make_pair(testdata[k].second,testdata[k].second*0.5)));
	}

	//H2H decrease
    cout<<"Decrease update..."<<endl;
	Graph g1=g;
	t1=std::chrono::high_resolution_clock::now();
	g1.H2HdecBat(testdataDec);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"H2H decrease "<<runT<<" "<<runT/testdata.size()<<endl;
    g1.CorCheckH2H(10000);

	//H2H increase
    cout<<"Increase update..."<<endl;
	Graph g2=g;
	t1=std::chrono::high_resolution_clock::now();
	g2.H2HincBatMT(testdataInc);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"H2H increase "<<runT<<" "<<runT/testdata.size()<<endl;
    g2.CorCheckH2H(10000);

	return 0;
}
