/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringGenerator.h"
#include "GraphClusteringTools.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Log.h"

namespace NetworKit {
Partition ClusteringGenerator::makeSingletonClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Partition ClusteringGenerator::makeOneClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.allToOnePartition();
	return zeta;
}

Partition ClusteringGenerator::makeRandomClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);

	zeta.setUpperBound(k);

	G.parallelForNodes([&](node v) {
		index c = Aux::Random::integer(k-1);
		zeta.addToSubset(c, v);
	});

	if (zeta.numberOfSubsets() != k) {
		WARN("random clustering does not contain k=",k," cluster: ",zeta.numberOfSubsets());
	}
	assert (GraphClusteringTools::isProperClustering(G, zeta));
	return zeta;
}

Partition ClusteringGenerator::makeContinuousBalancedClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound(); 
	Partition clustering(n);
	clustering.setUpperBound(k);

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
			clustering.addToSubset(block,v);
			++v;
		}
	}

	return clustering;
}

Partition ClusteringGenerator::makeNoncontinuousBalancedClustering(Graph &G, count k) {
	Partition clustering(G.upperNodeIdBound());
	clustering.setUpperBound(k);

	count i = 0;
	G.forNodes([&](node u) {
		clustering[u] = i % k;
		++i;
	});

	return clustering;
}

} /* namespace NetworKit */
