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

TEST(AdjacencyMatrixGTest, trySmallAdjacencyMatrix) {
	NetworKit::Graph graph(6);
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
	ASSERT_EQ(1, mat(0,0));
	ASSERT_EQ(1, mat(0,1));
	ASSERT_EQ(0, mat(0,2));
	ASSERT_EQ(0, mat(0,3));
	ASSERT_EQ(1, mat(0,4));
	ASSERT_EQ(0, mat(0,5));

	// third row
	ASSERT_EQ(0, mat(2,0));
	ASSERT_EQ(1, mat(2,1));
	ASSERT_EQ(0, mat(2,2));
	ASSERT_EQ(1, mat(2,3));
	ASSERT_EQ(0, mat(2,4));
	ASSERT_EQ(0, mat(2,5));

	// fifth row
	ASSERT_EQ(1, mat(4,0));
	ASSERT_EQ(1, mat(4,1));
	ASSERT_EQ(0, mat(4,2));
	ASSERT_EQ(1, mat(4,3));
	ASSERT_EQ(0, mat(4,4));
	ASSERT_EQ(0, mat(4,5));
}

TEST(AdjacencyMatrixGTest, tryAdjacencyMatrixOfLesmisGraph) {
	// read lesmis graph
	NetworKit::METISGraphReader graphReader;
	NetworKit::Graph graph = graphReader.read("input/lesmis.graph");

	// create AdjacencyMatrix
	AdjacencyMatrix mat(graph);

	mat.forElementsInRowOrder([&](const index row, const index column, const double value) {
		EXPECT_EQ(graph.weight(row, column), value);
	});
}


} /* namespace NetworKit */

