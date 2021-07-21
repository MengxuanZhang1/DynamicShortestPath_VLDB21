## Intro

This is the source code the VLDB 2020 paper "An Experimental Evaluation and Guideline for Path Finding in Weighted Dynamic Network". Please refer to the paper for the algorithm details.

## Algorithm  

The following codes contain the index construction and index update

1. CH.cpp: CH-P/CH-W
2. H2H.cpp: H2H
3. PLL.cpp: PLL
4. PathRetir.cpp: Path retrieval 
5. Founda.cpp: Index read/write, graph loading, efficiency evaluation, correctness check, DSP algorithms, index size, and query processing.
6. There are five main functions: mainCHW, mainCHP, mainH2H, mainPLL, mainDSP , which demonstrate the running example in the given network.

## Data

All the data are stored in the data folder.

One simple example network Graphsw10k5 is offered with node number of 1000, degree of 3 and the graph type is small-world. The first line stores the node number and edge number. The following lines store the nodeID1, nodeID2, length.

OD file stores all the queries, with the first line indicates the query number, and the following lines store the source and destination ID.

Order file stores all node orders, with the first line indicates the node number, and the following lines store the nodeID and node Order. The unreachable node's order is -1.

Update file stores all the updates, with the first line indicates the update number, and the following lines store the update edge's ID1, ID2, and the original weight. The new weights are generated in the main files according to the decrease / increase separately.

The constructed index files will also appear in the data folder.

## Dependency

g++ and boost.

All the codes are runnable after make.

