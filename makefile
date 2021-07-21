CXX=g++ -std=c++17
OPT=-O3

VLDB: CH.o Founda.o H2H.o PathRetri.o PLL.o mainCHP.o mainCHW.o mainDSP.o mainH2H.o mainPath.o mainPLL.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o CHP mainCHP.o CH.o Founda.o H2H.o PathRetri.o PLL.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o CHW mainCHW.o CH.o Founda.o H2H.o PathRetri.o PLL.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o DSP mainDSP.o CH.o Founda.o H2H.o PathRetri.o PLL.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o PLL mainPLL.o CH.o Founda.o H2H.o PathRetri.o PLL.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o H2H mainH2H.o CH.o Founda.o H2H.o PathRetri.o PLL.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o Path mainPath.o CH.o Founda.o H2H.o PathRetri.o PLL.o -lpthread -lboost_system -lboost_thread

mainCHP.o:mainCHP.cpp
	$(CXX) -g -c $(OPT) mainCHP.cpp
mainCHW.o:mainCHW.cpp
	$(CXX) -g -c $(OPT) mainCHW.cpp
mainDSP.o:mainDSP.cpp
	$(CXX) -g -c $(OPT) mainDSP.cpp
mainPLL.o:mainPLL.cpp
	$(CXX) -g -c $(OPT) mainPLL.cpp
mainH2H.o:mainH2H.cpp
	$(CXX) -g -c $(OPT) mainH2H.cpp
mainPath.o:mainPath.cpp
	$(CXX) -g -c $(OPT) mainPath.cpp
CH.o:CH.cpp
	$(CXX) -g -c $(OPT) CH.cpp
Founda.o:Founda.cpp
	$(CXX) -g -c $(OPT) Founda.cpp
H2H.o:H2H.cpp
	$(CXX) -g -c $(OPT) H2H.cpp
PathRetri.o:PathRetri.cpp
	$(CXX) -g -c $(OPT) PathRetri.cpp
PLL.o:PLL.cpp
	$(CXX) -g -c $(OPT) PLL.cpp

clean:
	rm*.o
	rm CHP
	rm CHW
	rm H2H
	rm PLL
	rm Path
	rm DSP
