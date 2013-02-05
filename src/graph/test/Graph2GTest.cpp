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
	int64_t n = 10;
	Graph2 G(n);

	EXPECT_EQ(n, G.numberOfNodes()) << "number of nodes should be " << n << " in the graph";


	// TODO: insert nodes and delete them, test again
}

TEST_F(Graph2GTest, testNumberOfEdges) {
	int64_t n = 10;
	Graph2 G(n);

	EXPECT_EQ(0, G.numberOfEdges()) << "no edges should be in the graph";

	G.insertEdge(1, 2);
	G.insertEdge(2, 3);
	G.insertEdge(1, 3);

	EXPECT_TRUE(G.hasEdge(1, 2));
	EXPECT_TRUE(G.hasEdge(2, 3));
	EXPECT_TRUE(G.hasEdge(1, 3));
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges should have been inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testIsEmpty) {
	Graph2 G(0);
	EXPECT_EQ(0, G.numberOfNodes()) << "no nodes should be in the graph";
	EXPECT_TRUE(G.isEmpty());

	// TODO: insert node, check if isEmpty is false, delete nodes again, check if true again
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
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges should have been inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
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

TEST_F(Graph2GTest, testConstNodeIteration) {
	// TODO:
}


TEST_F(Graph2GTest, testConstParallelNodeIteration) {
	// TODO:
}


TEST_F(Graph2GTest, testConstEdgeIteration) {
	// TODO:
}

TEST_F(Graph2GTest, testConstParallelEdgeIteration) {
	// TODO:
}

TEST_F(Graph2GTest, testConstNeighborIteration) {
	// TODO:
}

TEST_F(Graph2GTest, testConstParallelSumForNodes) {
	// TODO:
}

TEST_F(Graph2GTest, testConstNodePairIteration) {
	// TODO:
}


} /* namespace EnsembleClustering */
