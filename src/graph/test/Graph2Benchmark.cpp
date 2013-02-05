/*
 * Graph2Benchmark.cpp
 *
 *  Created on: 05.02.2013
 *      Author: cls
 */

#include "Graph2Benchmark.h"

namespace EnsembleClustering {

Graph2Benchmark::Graph2Benchmark() {
	// TODO Auto-generated constructor stub

}

Graph2Benchmark::~Graph2Benchmark() {
	// TODO Auto-generated destructor stub
}

TEST_F(Graph2Benchmark, graphConstruction) {

	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	count n = 1e+7;;
	Graph2 G(n);

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}


TEST_F(Graph2Benchmark, nodeIteration) {
	count n = 1e+7;;
	Graph2 G(n);


	std::vector<node> nodes(n, 0);

	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	G.forNodes([&](node v){
		nodes[v] = v;
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, parallelNodeIteration) {
	count n = 1e+7;;
	Graph2 G(n);


	std::vector<node> nodes(n, 0);

	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	G.parallelForNodes([&](node v){
		nodes[v] = v;
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}


TEST_F(Graph2Benchmark, nodePairIteration) {
	count n = 1e+4;;
	Graph2 G(n);


	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	count p = 0;
	G.forNodePairs([&](node u, node v){
		p += 1;
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

	EXPECT_EQ((n * (n-1)) / 2, p);
}


TEST_F(Graph2Benchmark, edgeInsertion) {

	count n = 1e+4;
	Graph2 G(n);

	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	G.forNodePairs([&](node u, node v){
		G.insertEdge(u, v);
	});


	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, edgeRemoval) {
	count n = 1e+4;
	Graph2 G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.insertEdge(u, v);
	});


	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	G.forNodePairs([&](node u, node v){
		G.removeEdge(u, v);
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, edgeIteration) {
	count n = 1e+4;
	Graph2 G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.insertEdge(u, v);
	});


	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	G.forEdges([&](node u, node v){

	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}

TEST_F(Graph2Benchmark, parallelEdgeIteration) {
	count n = 1e+4;
	Graph2 G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.insertEdge(u, v);
	});

	Aux::Timer run;
	INFO("[BEGIN]");
	run.start();

	// TODO:
	EXPECT_TRUE(false) << "TODO";

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}

TEST_F(Graph2Benchmark, parallelSumForNodes) {
	count n = 1e+7;
	Graph2 G(n);

	// TODO:
	EXPECT_TRUE(false) << "TODO";

}



TEST_F(Graph2Benchmark, nodeInsertion) {
	// TODO:
	EXPECT_TRUE(false) << "TODO";
}

TEST_F(Graph2Benchmark, nodeRemoval) {
	// TODO:
	EXPECT_TRUE(false) << "TODO";
}

} /* namespace EnsembleClustering */
