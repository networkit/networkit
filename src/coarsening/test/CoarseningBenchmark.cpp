/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "CoarseningBenchmark.h"

#include "../../generators/ErdosRenyiGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../coarsening/ClusterContracter.h"
#include "../../coarsening/PartitionCoarsening.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../auxiliary/Timer.h"

namespace NetworKit {

TEST_F(CoarseningBenchmark, benchmarkClusterContracter) {

	std::cout << "enter number of nodes: ";

	count n;
	std::cin >> n; 

	count redF = 100; // reduction factor
	DEBUG("generating graph with ", n, " nodes");
	auto gen = ErdosRenyiGenerator(n, 0.05);
	Graph G = gen.generate();

	DEBUG("generating partition");
	Clustering zeta(G.upperNodeIdBound());
	G.forNodes([&](node u){
		zeta[u] = u / redF;
	});

	DEBUG("number of subsets: ", zeta.numberOfClusters());

	INFO("sequential coarsening");
	Aux::Timer timer;
	ClusterContracter contracter;
	timer.start();
	auto result1 = contracter.run(G, zeta);
	timer.stop();
	INFO("sequential coarsening: ", timer.elapsedTag());

	Graph Gc1 = result1.first;

	INFO("parallel coarsening");
	PartitionCoarsening coarsening;
	timer.start();
	auto result2 = coarsening.run(G, zeta);
	timer.stop();
	Graph Gc2 = result2.first;
	INFO("parallel coarsening: ", timer.elapsedTag());


}





} /* namespace NetworKit */

#endif /*NOGTEST */
