/*
 * GraphGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef GRAPHGTEST_H_
#define GRAPHGTEST_H_

#include <gtest/gtest.h>

#include "../../aux/log.h"
#include "../Graph.h"
#include "../GraphGenerator.h"


namespace EnsembleClustering {

class GraphGTest: public testing::Test {

protected:

	GraphGenerator gen;

public:

	virtual void SetUp();

	virtual void TearDown();

};


TEST_F(GraphGTest, testEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	READ_ONLY_FORALL_EDGES_BEGIN(G) {
		node u = EDGE_SOURCE;
		node v = EDGE_DEST;
		if (u < v) {
			edgeCount += 1;
		}
	} READ_ONLY_FORALL_EDGES_END();

	EXPECT_EQ((n * (n-1)) / 2, edgeCount) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}


TEST_F(GraphGTest, testLambdaEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	G.forallEdges([&](node u, node v) {
		if (u < v) {
			edgeCount += 1;
		}
	},"readonly");

	EXPECT_EQ((n * (n-1)) / 2, edgeCount) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}



TEST_F(GraphGTest, testParallelLambdaEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	G.forallEdges([&](node u, node v) {
		if (u < v) {
			#pragma omp atomic update
			edgeCount += 1;
		}
	}, "parallel", "readonly");

	EXPECT_EQ((n * (n-1)) / 2, edgeCount) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}




TEST_F(GraphGTest, testLambdaEdgeModification) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	G.forallEdges([&](node u, node v){
		G.removeEdge(u, v);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been deleted";
}



TEST_F(GraphGTest, testLambdaNodeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t nodeCount = 0;
	G.forallNodes([&](node v) {
		nodeCount++;
	});

	EXPECT_EQ(n, nodeCount);
}


TEST_F(GraphGTest, testParallelLambdaNodeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t nodeCount = 0;
	G.forallNodes([&](node v) {
		#pragma omp atomic update
		nodeCount++;
	}, "parallel");

	EXPECT_EQ(n, nodeCount);
}


TEST_F(GraphGTest, testLambdaNeighborIteration) {

	Graph G;
	node v = 1;
	G.insertEdge(v, 2);
	G.insertEdge(v, 3);
	G.insertEdge(v, 4);

	int neighborCount = 0;
	G.forallNeighborsOf(v, [&](node w){
		neighborCount += 1;
	});

	EXPECT_EQ(3, neighborCount) << "node v has 3 neighbors";
}


TEST_F(GraphGTest, testLambdaIncidentEdgeIteration) {

	Graph G;
	node v = 1;
	G.insertEdge(v, 2);
	G.insertEdge(v, 3);
	G.insertEdge(v, 4);

	int edgeCount = 0;
	G.forallEdgesOf(v, [&](node v, node w){
		edgeCount += 1;
	});

	EXPECT_EQ(3, edgeCount) << "node v has 3 incident edges";
}


} /* namespace EnsembleClustering */
#endif /* GRAPHGTEST_H_ */
