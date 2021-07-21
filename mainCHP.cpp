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
	Graph g;
	string destiGraph=DesFile+"Graphsw10k5";
	string indexfiledec=DesFile+"CHP";
	string indexfileinc=DesFile+"CHPinc";
	string orderfile=DesFile+"Order";
	//read index update data
	string updateFile=DesFile+"Update";
	vector<pair<pair<int,int>,int>> testdata;
	g.ReadIns(updateFile, testdata);

	//*********************for decrease**************************//
	g.ReadGraph(destiGraph);
	//index construction
	t1=std::chrono::high_resolution_clock::now();
	g.StartCHPOrderMT(indexfiledec,orderfile);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CHP decrease construction time "<<runT<<endl;
	//index size & efficiency
	string ODfile=DesFile+"OD";
	g.CHPindexsize();
	g.EffiCheckCHPorder(ODfile);
	g.CorCheckCHPorder();
	//index update
	Graph g1=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g1.CHPdec(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*0.5);
		//g.CorCheckCHPorder();
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CHP decrease "<<runT<<" "<<runT/testdata.size()<<endl;

	//*********************for increase**************************//
	g.ReadGraphForCHPinc(destiGraph);
	//index construction
	t1=std::chrono::high_resolution_clock::now();
	g.StartCHPOrderMTWP(indexfileinc,orderfile);
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CHP increase construction time "<<runT<<endl;
	//index size & efficiency
	g.CHPindexsize();
	g.CHPWitPathsize();
	g.EffiCheckCHPorder(ODfile);
	g.CorCheckCHPorder();
	//index update
	Graph g2=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g2.CHPinc(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*10);
		//g.CorCheckCHPorder();
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"CHP increase "<<runT<<" "<<runT/testdata.size()<<endl;

	return 0;
}
