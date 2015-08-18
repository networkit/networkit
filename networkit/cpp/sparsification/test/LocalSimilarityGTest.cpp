/*
 * LocalSimilarityGTest.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "LocalSimilarityGTest.h"

#include "../LocalSimilarityAttributizer.h"
#include "../../edgescores/ChibaNishizekiTriangleCounter.h"


namespace NetworKit {

TEST_F(LocalSimilarityGTest, testAttributeSimple) {
	Graph g(4);

	g.addEdge(0, 1);
	g.addEdge(0, 3);
	g.addEdge(0, 2);
	g.addEdge(1, 2);
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleCounter(g);
	std::vector<count> triangles = triangleCounter.getAttribute();
	LocalSimilarityAttributizer localSim(g, triangles);
	std::vector<double> exp = localSim.getAttribute();

	EXPECT_DOUBLE_EQ(1.0, exp[g.edgeId(0, 1)]);
	EXPECT_NEAR(0.36907025, exp[g.edgeId(0, 2)], 1e-7);
	EXPECT_DOUBLE_EQ(1.0, exp[g.edgeId(0, 3)]);
	EXPECT_DOUBLE_EQ(1.0, exp[g.edgeId(1, 2)]);
}

}
/* namespace NetworKit */

#endif /*NOGTEST */
