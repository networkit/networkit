/*
 * DistMeasureTest.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>
#include <cstdio>

#include "../../graph/Graph.h"
#include "../../viz/PostscriptWriter.h"
#include "../../io/METISGraphReader.h"
#include "../../io/DibapGraphReader.h"
#include "../../structures/Partition.h"
#include "../../community/Modularity.h"
#include "../AlgebraicDistanceIndex.h"

#include "../../auxiliary/Log.h"

namespace NetworKit {

class DistMeasuresGTest: public testing::Test {};

TEST_F(DistMeasuresGTest, testAlgebraicDistanceIndex) {
	Graph G(42);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v);
	});

	count numSystems = 2;
	count numIterations = 200;
	double omega = 0.5;
	AlgebraicDistanceIndex ad(G, numSystems, numIterations, omega);
	ad.preprocess();

	double adSum = 0.0;
	G.forNodePairs([&](node u, node v){
		adSum += ad.run(u, v);
	});

	DEBUG("sum of algebraic distances: " , adSum);
	EXPECT_GE(1e-12, adSum) << "algebraic distances should converge towards zero";
}

} /* namespace NetworKit */

