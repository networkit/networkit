/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "CoarseningBenchmark.h"

#include "../../generators/ErdosRenyiGenerator.h"
#include "../../community/ClusteringGenerator.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../community/ClusteringGenerator.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Log.h"

namespace NetworKit {

TEST_F(CoarseningBenchmark, benchmarkCoarsening) {
	std::cout << "enter number of nodes: ";
	count n;
	std::cin >> n;
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
	ParallelPartitionCoarsening coarsening(false);
	timer.start();
	auto result2 = coarsening.run(G, zeta);
	timer.stop();
	Graph Gc2 = result2.first;
	INFO("parallel coarsening: ", timer.elapsedTag());
	EXPECT_EQ(k, Gc2.numberOfNodes());

	INFO("parallel coarsening using GraphBuilder");
	ParallelPartitionCoarsening gbCoarsening(true);
	timer.start();
	auto result3 = gbCoarsening.run(G, zeta);
	timer.stop();
	Graph Gc3 = result3.first;
	INFO("parallel coarsening: ", timer.elapsedTag());
	EXPECT_EQ(k, Gc3.numberOfNodes());

}



} /* namespace NetworKit */

#endif /*NOGTEST */
