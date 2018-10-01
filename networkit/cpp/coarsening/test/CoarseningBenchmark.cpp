/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>

#include "../../generators/ErdosRenyiGenerator.h"
#include "../../community/ClusteringGenerator.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../community/ClusteringGenerator.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Log.h"

namespace NetworKit {

class CoarseningBenchmark: public testing::Test {};

TEST_F(CoarseningBenchmark, benchmarkCoarsening) {
	count n = 10000;
	count redF = 100; // reduction factor
	count k = n/redF;
	DEBUG("generating graph with ", n, " nodes");
	auto gen = ErdosRenyiGenerator(n, 0.05);
	Graph G = gen.generate();

	DEBUG("generating random partition");
	ClusteringGenerator clusteringGen;
	Partition zeta = clusteringGen.makeRandomClustering(G, k);

	//count k = zeta.numberOfSubsets();
	DEBUG("number of subsets: ", k);

	Aux::Timer timer;
	INFO("parallel coarsening");
	timer.start();
	ParallelPartitionCoarsening coarsening(G, zeta, false);
	coarsening.run();
	timer.stop();
	Graph Gc2 = coarsening.getCoarseGraph();
	INFO("parallel coarsening: ", timer.elapsedTag());
	EXPECT_EQ(k, Gc2.numberOfNodes());

	INFO("parallel coarsening using GraphBuilder");
	timer.start();
	ParallelPartitionCoarsening gbCoarsening(G, zeta, true);
	gbCoarsening.run();
	timer.stop();
	Graph Gc3 = gbCoarsening.getCoarseGraph();
	INFO("parallel coarsening: ", timer.elapsedTag());
	EXPECT_EQ(k, Gc3.numberOfNodes());
}

} /* namespace NetworKit */

