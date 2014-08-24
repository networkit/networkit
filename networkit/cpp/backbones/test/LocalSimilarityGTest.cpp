/*
 * LocalSimilarityGTest.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "LocalSimilarityGTest.h"

#include "../Backbones.h"
#include "../LocalSimilarityAttributizer.h"


namespace NetworKit {

TEST_F(LocalSimilarityGTest, testSimilarityCalculation) {
	Graph g(10);

	g.addEdge(0, 2);
	g.addEdge(0, 3);
	g.addEdge(0, 4);
	g.addEdge(0, 5);
	g.addEdge(0, 9);
	g.addEdge(1, 4);
	g.addEdge(1, 5);
	g.addEdge(1, 6);
	g.addEdge(1, 7);
	g.addEdge(1, 8);
	g.addEdge(1, 9);
	g.indexEdges();

	LocalSimilarityAttributizer localSim;

	//Common neighbors of 0 and 1: 4, 5, 9 --> jaccard measure of 3/8.
	EXPECT_DOUBLE_EQ(0.375, localSim.getSimilarity(g, 0, 1));
	EXPECT_DOUBLE_EQ(0.375, localSim.getSimilarity(g, 1, 0));
}

TEST_F(LocalSimilarityGTest, testAttributeSimple) {
	Graph g(4);

	g.addEdge(0, 1);
	g.addEdge(0, 3);
	g.addEdge(0, 2);
	g.addEdge(1, 2);
	g.indexEdges();

	LocalSimilarityAttributizer localSim;
	std::vector<double> exp = localSim.getAttribute(g, std::vector<int>(g.upperEdgeIdBound()));

	EXPECT_DOUBLE_EQ(0.0, exp[g.edgeId(0, 1)]);
	EXPECT_NEAR(0.63092975, exp[g.edgeId(0, 2)], 1e-7);
	EXPECT_DOUBLE_EQ(0.0, exp[g.edgeId(0, 3)]);
	EXPECT_DOUBLE_EQ(0.0, exp[g.edgeId(1, 2)]);
}

}
/* namespace NetworKit */

#endif /*NOGTEST */
