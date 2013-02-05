/*
 * Graph2GTest.cpp
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#include "Graph2GTest.h"

namespace EnsembleClustering {

Graph2GTest::Graph2GTest() {
	// TODO Auto-generated constructor stub

}

Graph2GTest::~Graph2GTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(Graph2GTest, testNumberOfNodes) {
	// TODO:
}

TEST_F(Graph2GTest, testNumberOfEdges) {
	// TODO:
}

TEST_F(Graph2GTest, testIsEmpty) {
	// TODO:
}

TEST_F(Graph2GTest, testEdgeInsertionAndRemoval) {

	int64_t n = 10;
	Graph2 G(n);

	EXPECT_EQ(0, G.numberOfEdges()) << "no edges should be in the graph";

	G.insertEdge(1, 2);
	G.insertEdge(2, 3);
	G.insertEdge(1, 3);

	EXPECT_TRUE(G.hasEdge(1, 2));
	EXPECT_TRUE(G.hasEdge(2, 3));
	EXPECT_TRUE(G.hasEdge(1, 3));
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges were inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges were removed";
}

TEST_F(Graph2GTest, testNodeIteration) {
	// TODO:
}


TEST_F(Graph2GTest, testParallelNodeIteration) {
	// TODO:
}


TEST_F(Graph2GTest, testEdgeIteration) {
	// TODO:
}

TEST_F(Graph2GTest, testParallelEdgeIteration) {
	// TODO:
}

TEST_F(Graph2GTest, testNeighborIteration) {
	// TODO:
}

TEST_F(Graph2GTest, testParallelSumForNodes) {
	// TODO:
}

TEST_F(Graph2GTest, testNodePairIteration) {
	// TODO:
}

} /* namespace EnsembleClustering */
