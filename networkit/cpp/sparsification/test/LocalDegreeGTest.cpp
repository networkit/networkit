/*
 * LocalDegreeGTest.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "LocalDegreeGTest.h"
#include "../LocalDegreeAttributizer.h"
#include <math.h>

namespace NetworKit {

TEST_F(LocalDegreeGTest, testAttributeSimple) {
	Graph g(22);
	g.addEdge(0, 1);
	g.addEdge(0, 2);
	g.addEdge(2, 3);
	g.addEdge(3, 4);
	g.addEdge(2, 4);
	g.addEdge(2, 5);
	g.addEdge(2, 6);
	g.addEdge(2, 7);
	g.addEdge(4, 7);
	g.addEdge(5, 8);
	g.addEdge(5, 9);
	g.addEdge(5, 10);
	g.addEdge(5, 11);
	g.addEdge(5, 12);
	g.addEdge(6, 13);
	g.addEdge(6, 14);
	g.addEdge(6, 15);
	g.addEdge(6, 16);
	g.addEdge(7, 17);
	g.addEdge(7, 18);
	g.addEdge(7, 19);
	g.addEdge(3, 20);
	g.addEdge(3, 21);
	g.indexEdges();

	LocalDegreeAttributizer localDegree(g);
	std::vector<double> attribute = localDegree.getAttribute();

	EXPECT_DOUBLE_EQ(LocalDegreeGTest::getScore(g, 0, 1, 1, 2), attribute[g.edgeId(0, 1)]);
	EXPECT_DOUBLE_EQ(LocalDegreeGTest::getScore(g, 2, 4, 1, 4), attribute[g.edgeId(2, 4)]);
	EXPECT_DOUBLE_EQ(LocalDegreeGTest::getScore(g, 4, 7, 2, 2), attribute[g.edgeId(4, 7)]);
}

/***
Calculates the LD score for an edge.
@param g the graph
@param x first node
@param y second node
@param rankX rank of y in the neighborhood of y (1-based)
@param rankY rank of y in the neighborhood of x (1-based)
**/
double LocalDegreeGTest::getScore(const Graph& g, node x, node y, count rankX, count rankY) {
	//Special case: degree one
	if (g.degree(x) == 1 || g.degree(y) == 1)
		return 1;

	//Use only the edge that leads to the higher-degree-node!
	if (g.degree(x) > g.degree(y))
		return 1 - log(rankX) / log(g.degree(y));

	if (g.degree(x) < g.degree(y))
		return 1 - log(rankY) / log(g.degree(x));

	return -1;
}

}
/* namespace NetworKit */

#endif /*NOGTEST */
