/*
 * AdjacencyMatrixGTest.cpp
 *
 *  Created on: 02.04.2014
 *      Author: Michael
 */

#include "AdjacencyMatrixGTest.h"

namespace NetworKit {

AdjacencyMatrixGTest::AdjacencyMatrixGTest() {
}

AdjacencyMatrixGTest::~AdjacencyMatrixGTest() {
}

TEST(AdjacencyMatrixGTest, testSmallAdjacencyMatrix) {
	Graph graph(6);
	graph.addEdge(0,0);
	graph.addEdge(0,1);
	graph.addEdge(0,4);
	graph.addEdge(1,2);
	graph.addEdge(1,4);
	graph.addEdge(2,3);
	graph.addEdge(3,4);
	graph.addEdge(3,5);

	AdjacencyMatrix mat(graph);

	// first row
	EXPECT_EQ(1, mat(0,0));
	EXPECT_EQ(1, mat(0,1));
	EXPECT_EQ(0, mat(0,2));
	EXPECT_EQ(0, mat(0,3));
	EXPECT_EQ(1, mat(0,4));
	EXPECT_EQ(0, mat(0,5));

	// third row
	EXPECT_EQ(0, mat(2,0));
	EXPECT_EQ(1, mat(2,1));
	EXPECT_EQ(0, mat(2,2));
	EXPECT_EQ(1, mat(2,3));
	EXPECT_EQ(0, mat(2,4));
	EXPECT_EQ(0, mat(2,5));

	// fifth row
	EXPECT_EQ(1, mat(4,0));
	EXPECT_EQ(1, mat(4,1));
	EXPECT_EQ(0, mat(4,2));
	EXPECT_EQ(1, mat(4,3));
	EXPECT_EQ(0, mat(4,4));
	EXPECT_EQ(0, mat(4,5));


	// directed, weighted graph
	Graph dGraph(4, true, true);
	dGraph.addEdge(0,1,2);
	dGraph.addEdge(0,0, 42);
	dGraph.addEdge(2,3,-3);
	dGraph.addEdge(3,2,5);

	mat = AdjacencyMatrix(dGraph);
	ASSERT_EQ(dGraph.numberOfNodes(), mat.numberOfRows());
	ASSERT_EQ(dGraph.numberOfNodes(), mat.numberOfColumns());

	EXPECT_EQ(2, mat(0,1));
	EXPECT_EQ(0, mat(1,0));
	EXPECT_EQ(42, mat(0,0));
	EXPECT_EQ(-3, mat(2,3));
	EXPECT_EQ(5, mat(3,2));
}

TEST(AdjacencyMatrixGTest, testAdjacencyMatrixOfLesmisGraph) {
	// read lesmis graph
	METISGraphReader graphReader;
	Graph graph = graphReader.read("input/lesmis.graph");

	// create AdjacencyMatrix
	AdjacencyMatrix mat(graph);

	mat.forNonZeroElementsInRowOrder([&](const index row, const index column, const double value) {
		EXPECT_EQ(graph.weight(row, column), value);
	});
}


} /* namespace NetworKit */

