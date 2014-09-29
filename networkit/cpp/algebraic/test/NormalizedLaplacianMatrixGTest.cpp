/*
 * NormalizedLaplacianMatrixGTest.cpp
 *
 *  Created on: 25.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "NormalizedLaplacianMatrixGTest.h"


namespace NetworKit {

NormalizedLaplacianMatrixGTest::NormalizedLaplacianMatrixGTest() {
}

NormalizedLaplacianMatrixGTest::~NormalizedLaplacianMatrixGTest() {
}

TEST(NormalizedLaplacianMatrixGTest, testSmallNormalizedLaplacianMatrix) {
	NetworKit::Graph graph(7);
	graph.addEdge(0, 1);
	graph.addEdge(0, 4);
	graph.addEdge(1, 4);
	graph.addEdge(1, 2);
	graph.addEdge(2, 3);
	graph.addEdge(3, 4);
	graph.addEdge(3, 5);

	NormalizedLaplacianMatrix normalizedLaplacianMatrix(graph);
	ASSERT_EQ(graph.numberOfNodes(), normalizedLaplacianMatrix.numberOfRows());
	ASSERT_EQ(graph.numberOfNodes(), normalizedLaplacianMatrix.numberOfColumns());

	EXPECT_EQ(1, normalizedLaplacianMatrix(0,0));
	EXPECT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(0,1));
	EXPECT_EQ(0, normalizedLaplacianMatrix(0,2));
	EXPECT_EQ(0, normalizedLaplacianMatrix(0,3));
	EXPECT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(0,4));
	EXPECT_EQ(0, normalizedLaplacianMatrix(0,5));
	EXPECT_EQ(0, normalizedLaplacianMatrix(0,6));
	EXPECT_EQ(1, normalizedLaplacianMatrix(1,1));
	EXPECT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(1,2));
	EXPECT_EQ(0, normalizedLaplacianMatrix(1,3));
	EXPECT_EQ(-1.0 / 3.0, normalizedLaplacianMatrix(1,4));
	EXPECT_EQ(0, normalizedLaplacianMatrix(1,5));
	EXPECT_EQ(0, normalizedLaplacianMatrix(1,6));
	EXPECT_EQ(1, normalizedLaplacianMatrix(2,2));
	EXPECT_EQ(-1.0 / sqrt(2.0 * 3.0), normalizedLaplacianMatrix(2,3));
	EXPECT_EQ(0, normalizedLaplacianMatrix(2,4));
	EXPECT_EQ(0, normalizedLaplacianMatrix(2,5));
	EXPECT_EQ(0, normalizedLaplacianMatrix(2,6));
	EXPECT_EQ(1, normalizedLaplacianMatrix(3,3));
	EXPECT_EQ(-1.0 / 3.0, normalizedLaplacianMatrix(3,4));
	EXPECT_EQ(-1.0 / sqrt(3.0), normalizedLaplacianMatrix(3,5));
	EXPECT_EQ(0, normalizedLaplacianMatrix(3,6));
	EXPECT_EQ(1, normalizedLaplacianMatrix(4,4));
	EXPECT_EQ(0, normalizedLaplacianMatrix(4,5));
	EXPECT_EQ(0, normalizedLaplacianMatrix(4,6));
	EXPECT_EQ(1, normalizedLaplacianMatrix(5,5));
	EXPECT_EQ(0, normalizedLaplacianMatrix(5,6));
	EXPECT_EQ(0, normalizedLaplacianMatrix(6,6));
}

TEST(NormalizedLaplacianMatrixGTest, testNormalizedLaplacianMatrixOfLesmisGraph) {
	// read lesmis graph
	NetworKit::METISGraphReader graphReader;
	NetworKit::Graph graph = graphReader.read("input/lesmis.graph");

	// create NormalizedLaplacianMatrix
	NormalizedLaplacianMatrix mat(graph);

	mat.forElementsInRowOrder([&](const index row, const index column, const double value){
		if (row == column) {
			if (graph.weightedDegree(row) != 0) {
				if (graph.isWeighted()) {
					EXPECT_EQ(1-graph.weight(row, row)/graph.weightedDegree(row), value);
				} else {
					EXPECT_EQ(1, value);
				}
			}
		} else {
			EXPECT_EQ(-graph.weight(row, column)/sqrt(graph.weightedDegree(row)*graph.weightedDegree(column)), value);
		}
	});
}

} /* namespace NetworKit */

