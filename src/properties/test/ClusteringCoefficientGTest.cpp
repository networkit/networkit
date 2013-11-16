/*
 * ClusteringCoefficientGTest.cpp
 *
 *  Created on: Nov 16, 2013
 *      Author: lbarth, dweiss
 */
#ifndef NOGTEST

#include "ClusteringCoefficientGTest.h"

namespace NetworKit {

ClusteringCoefficientGTest::ClusteringCoefficientGTest() {
	// TODO Auto-generated constructor stub
}

ClusteringCoefficientGTest::~ClusteringCoefficientGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(ClusteringCoefficientGTest, testRunByCall) {
	// construct graph
	Graph g;
	for (count i = 0; i < 20; i++) {
		g.addNode();
	}
	g.addEdge(0,1,0);
	g.addEdge(1,2,0);
	g.addEdge(2,4,0);
	g.addEdge(4,8,0);
	g.addEdge(8,16,0);
	g.addEdge(16,19,0);

	g.addEdge(3,5,0);
	g.addEdge(5,6,0);
	g.addEdge(6,7,0);
	g.addEdge(7,9,0);

	g.addEdge(10,11,0);
	g.addEdge(10,18,0);
	g.addEdge(10,12,0);
	g.addEdge(18,17,0);

	g.addEdge(13,14,0);

	// initialize ClusteringCoefficient
	ClusteringCoefficient cc;

	// compute avg. local CC
	double avgLocalCC = cc.avgLocal(g);

	EXPECT_TRUE(avgLocalCC == 2.0);
}


} /* namespace NetworKit */

#endif /*NOGTEST */

