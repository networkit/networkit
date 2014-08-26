/*
 * ExperimentalGTest.cpp
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "ExperimentalGTest.h"

#include "../TopDegreeAttributizer.h"


namespace NetworKit {

TEST_F(ExperimentalGTest, testSimilarityCalculation) {
	Graph g(8);

	g.addEdge(0,2);
	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(3,4);
	g.addEdge(4,5);
	g.addEdge(5,6);
	g.addEdge(4,6);
	g.addEdge(4,7);
	g.indexEdges();

	TopDegreeAttributizer tdAttributizer;
	std::vector<count> expected = tdAttributizer.getAttribute(g, std::vector<int>(0));
	EXPECT_EQ(1, expected[g.edgeId(0,2)]);
	EXPECT_EQ(1, expected[g.edgeId(1,2)]);
	EXPECT_EQ(1, expected[g.edgeId(2,3)]);
	EXPECT_EQ(1, expected[g.edgeId(3,4)]);
	EXPECT_EQ(1, expected[g.edgeId(4,5)]);
	EXPECT_EQ(2, expected[g.edgeId(5,6)]);
	EXPECT_EQ(1, expected[g.edgeId(4,6)]);
	EXPECT_EQ(1, expected[g.edgeId(4,7)]);
}

}

/* namespace NetworKit */

#endif /*NOGTEST */
