/*
 * ClusteringGenerator2.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringGenerator2.h"
#include "../auxiliary/Random.h"
#include "../structures/Partition.h"

namespace NetworKit {

ClusteringGenerator2::ClusteringGenerator2() {

}

ClusteringGenerator2::~ClusteringGenerator2() {

}

Partition ClusteringGenerator2::makeSingletonClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Partition ClusteringGenerator2::makeOneClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.toSingleton(0);
	index one = zeta[0];
	G.forNodes([&](node v){
		zeta.addToSubset(one, v);
	});
	return zeta;
}

Partition ClusteringGenerator2::makeRandomClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);

	/*for (uint64_t i = 0; i < k; ++i) {
		zeta.toSingleton(i);
	}*/
	zeta.setUpperBound(k);

	G.parallelForNodes([&](node v){ //parallelForNodes
		index c = Aux::Random::integer(k-1); 
		zeta.addToSubset(c, v);
	});

	//assert (zeta.isProper(G)); #not existent in partition
	return zeta;
}

Partition ClusteringGenerator2::makeContinuousBalancedClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
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
		clustering.toSingleton(v);
		index p = clustering[v];
		++v;
		while (v < blockSize[block]) {
			clustering.addToSubset(p,v);
			++v;
		}
	}

	return clustering;
}

} /* namespace NetworKit */
