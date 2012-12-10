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
	Clustering zeta(G.numberOfNodes());
	for (node v = G.firstNode(); v <= G.numberOfNodes(); ++v) {
		zeta.toSingleton(v);
	}
	return zeta;
}

Clustering& ClusteringGenerator::makeOneClustering(const Graph& G) {
	Clustering zeta(G.numberOfNodes());
	cluster one = 1;
	for (node v = G.firstNode(); v <= G.numberOfNodes(); ++v) {
		zeta.addToCluster(one, v);
	}
	return zeta;
}

} /* namespace EnsembleClustering */
