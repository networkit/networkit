/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#include "ClusteringGenerator.h"

namespace EnsembleClustering {

ClusteringGenerator::ClusteringGenerator() {
	// TODO Auto-generated constructor stub

}

ClusteringGenerator::~ClusteringGenerator() {
	// TODO Auto-generated destructor stub
}

Clustering& ClusteringGenerator::makeSingletonClustering(const Graph& G) {
	int64_t n = G.numberOfNodes();
	Clustering* zeta = new Clustering(n);
	zeta->allToSingletons();
	return *zeta;
}

Clustering& ClusteringGenerator::makeOneClustering(const Graph& G) {
	int64_t n = G.numberOfNodes();
	Clustering* zeta = new Clustering(n);
	cluster one = 1;
	for (node v = G.firstNode(); v <= n; ++v) {
		zeta->addToCluster(one, v);
	}
	return *zeta;
}

} /* namespace EnsembleClustering */
