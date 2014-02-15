/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringGenerator.h"
#include "../auxiliary/Random.h"
#include "GraphClusteringTools.h"

namespace NetworKit {

ClusteringGenerator::ClusteringGenerator() {

}

ClusteringGenerator::~ClusteringGenerator() {

}

Partition ClusteringGenerator::makeSingletonClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Partition ClusteringGenerator::makeOneClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.setUpperBound(1);
	G.forNodes([&](node v){
		zeta.addToSubset(0, v); //TODO not very nice...
	});
	return zeta;
}

Partition ClusteringGenerator::makeRandomClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);

	zeta.setUpperBound(k-1);

	G.parallelForNodes([&](node v){
		index c = Aux::Random::integer(k-1);
		zeta.addToSubset(c, v);
	});

	assert (GraphClusteringTools::isProperClustering(G, zeta));
	return zeta;
}

Partition ClusteringGenerator::makeContinuousBalancedClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound(); // FIXME: upper Node ID bound is actually not the right way to do this.
	Partition clustering(n);

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
			clustering.addToSubset(block,v);//clustering[v] = block;
			++v;
		}
	}

	return clustering;
}

} /* namespace NetworKit */
