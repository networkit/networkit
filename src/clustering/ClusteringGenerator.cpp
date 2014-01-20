/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringGenerator.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

ClusteringGenerator::ClusteringGenerator() {

}

ClusteringGenerator::~ClusteringGenerator() {

}

Clustering ClusteringGenerator::makeSingletonClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Clustering zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Clustering ClusteringGenerator::makeOneClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Clustering zeta(n);
	cluster one = zeta.addCluster();
	G.forNodes([&](node v){
		zeta.addToCluster(one, v);
	});
	return zeta;
}

Clustering ClusteringGenerator::makeRandomClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
	Clustering zeta(n);

	for (uint64_t i = 0; i < k; ++i) {
		zeta.addCluster();
	}

	G.parallelForNodes([&](node v){
		cluster c = Aux::Random::integer(k-1);
		zeta.addToCluster(c, v);
	});

	assert (zeta.isProper(G));
	return zeta;
}

Clustering ClusteringGenerator::makeContinuousBalancedClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
	Clustering clustering(n);

	std::vector<count> blockSize(k, 0);

	// compute block sizes
	for (index block = 0; block < k; ++block) {
		blockSize[block] = n / k + (n % k > block);
	}

	// compute prefix sums of block sizes
	for (index block = 1; block < k; ++block) {
		blockSize[block] += blockSize[block-1];
	}

	// fill clustering according to blocks
	node v = 0;
	for (index block = 0; block < k; ++block) {
		while (v < blockSize[block]) {
			clustering[v] = block;
			++v;
		}
	}

	return clustering;
}

} /* namespace NetworKit */
