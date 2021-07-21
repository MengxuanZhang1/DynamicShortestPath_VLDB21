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
	string PLLfile=DesFile+"PLL";
	string PLLpointfile=DesFile+"PLLPoint";
	g.ReadGraph(destiGraph);

	//index construction
	g.StartPLL(PLLfile, PLLpointfile, orderfile);
	g.CorCheckPLL();

	//index contruction without write & read
	/*t1=std::chrono::high_resolution_clock::now();
	g.PLLcon(orderfile);//not read & write index
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL construction time "<<runT<<endl;*/

	//index size & efficiency
	string ODfile=DesFile+"OD";
	g.PLLindexsize();
	g.PrunePointsize();
	g.EffiCheckPLL(ODfile);

	//index update
	string updateFile=DesFile+"Update";
	vector<pair<pair<int,int>,int>> testdata;
	g.ReadIns(updateFile, testdata);

	//PLL-p decrease
	Graph g1=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g1.PSLdec(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*0.5);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL-p decrease "<<runT<<" "<<runT/testdata.size()<<endl;

	//PLL-p increase
	Graph g2=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g2.PSLinc(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*1.5);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL-p increase "<<runT<<" "<<runT/testdata.size()<<endl;

	//PLL-s decrease
	Graph g3=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g3.PLLdec(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*0.5);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL-s decrease "<<runT<<" "<<runT/testdata.size()<<endl;

	//PLL-s increase
	Graph g4=g;
	t1=std::chrono::high_resolution_clock::now();
	for(int k=0;k<testdata.size();k++){
		g4.PLLinc(testdata[k].first.first, testdata[k].first.second, testdata[k].second, testdata[k].second*1.5);
	}
	t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"PLL-s increase "<<runT<<" "<<runT/testdata.size()<<endl;


	return 0;
}



