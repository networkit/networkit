/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringGenerator.h"

namespace NetworKit {

ClusteringGenerator::ClusteringGenerator() {
	// TODO Auto-generated constructor stub

}

ClusteringGenerator::~ClusteringGenerator() {
	// TODO Auto-generated destructor stub
}

Clustering ClusteringGenerator::makeSingletonClustering(Graph& G) {
	count n = G.numberOfNodes();
	Clustering zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Clustering ClusteringGenerator::makeOneClustering(Graph& G) {
	count n = G.numberOfNodes();
	Clustering zeta(n);
	cluster one = zeta.addCluster();
	G.forNodes([&](node v){
		zeta.addToCluster(one, v);
	});
	return zeta;
}

Clustering ClusteringGenerator::makeRandomClustering(Graph& G, count k) {
	// random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, k-1);
	//
	count n = G.numberOfNodes();
	Clustering zeta(n);

	for (int64_t i = 0; i < k; ++i) {
		zeta.addCluster();
	}

	G.parallelForNodes([&](node v){
		cluster c = dis(gen);
		zeta.addToCluster(c, v);
	});

	assert (zeta.isProper(G));
	return zeta;
}

Clustering ClusteringGenerator::makeContinuousBalancedClustering(Graph& G, count k) {
	count n = G.numberOfNodes();
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
