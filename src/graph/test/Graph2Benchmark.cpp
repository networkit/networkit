/*
 * Graph2Benchmark.cpp
 *
 *  Created on: 05.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "Graph2Benchmark.h"

namespace NetworKit {

Graph2Benchmark::Graph2Benchmark() {
	// TODO Auto-generated constructor stub

}

Graph2Benchmark::~Graph2Benchmark() {
	// TODO Auto-generated destructor stub
}

TEST_F(Graph2Benchmark, graphConstruction) {
	count n = 1e+7;;

	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	Graph G(n);

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}


TEST_F(Graph2Benchmark, nodeIteration) {
	count n = 1e+7;;
	Graph G(n);


	std::vector<node> nodes(n, 0);

	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	G.forNodes([&](node v){
		nodes[v] = v;
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, parallelNodeIteration) {
	count n = 1e+7;;
	Graph G(n);


	std::vector<node> nodes(n, 0);

	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	G.parallelForNodes([&](node v){
		nodes[v] = v;
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}


TEST_F(Graph2Benchmark, nodePairIteration) {
	count n = 1e+4;;
	Graph G(n);


	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
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
	Graph G(n);

	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});


	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, parallelEdgeInsertion) {

	count n = 1e+4;
	Graph G(n);

	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	// TODO: test parallel edge insertion
	EXPECT_TRUE(false) << "TODO";

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, edgeRemoval) {
	count n = 1e+4;
	Graph G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});


	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	G.forNodePairs([&](node u, node v){
		G.removeEdge(u, v);
	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, parallelEdgeRemoval) {
	count n = 1e+4;
	Graph G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});


	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	// TODO: test parallel edge removal
	EXPECT_TRUE(false) << "TODO";

	run.stop();
	INFO("[DONE]" << run.elapsedTag());

}

TEST_F(Graph2Benchmark, edgeIteration) {
	count n = 1e+4;
	Graph G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});


	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	G.forEdges([&](node u, node v){

	});

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}

TEST_F(Graph2Benchmark, parallelEdgeIteration) {
	count n = 1e+4;
	Graph G(n);

	// insert edges
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});

	Aux::Timer run;
	INFO("[BEGIN] (n=" << n << ")");
	run.start();

	// TODO: benchmark parallel edge iteration
	EXPECT_TRUE(false) << "TODO";

	run.stop();
	INFO("[DONE]" << run.elapsedTag());
}

TEST_F(Graph2Benchmark, parallelSumForNodes) {
	count n = 1e+7;
	Graph G(n);

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

} /* namespace NetworKit */

#endif /*NOGTEST */
