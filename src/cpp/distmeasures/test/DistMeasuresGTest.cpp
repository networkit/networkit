/*
 * DistMeasureTest.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "DistMeasuresGTest.h"
#include "../../auxiliary/Log.h"

namespace NetworKit {

DistMeasuresGTest::DistMeasuresGTest() {

}

DistMeasuresGTest::~DistMeasuresGTest() {

}

TEST_F(DistMeasuresGTest, testAlgebraicDistance) {
	GraphGenerator gen;
	Graph G = gen.makeCompleteGraph(42);

	count numSystems = 2;
	count numIterations = 200;
	double omega = 0.5;
	NodeDistance* ad = new AlgebraicDistance(G, numSystems, numIterations, omega);
	ad->preprocess();

	double adSum = 0.0;
	G.forNodePairs([&](node u, node v){
		adSum += ad->distance(u, v);
	});

	INFO("sum of algebraic distances: " , adSum);
	EXPECT_GE(1e-12, adSum) << "algebraic distances should converge towards zero";

}


} /* namespace NetworKit */

#endif /*NOGTEST */
