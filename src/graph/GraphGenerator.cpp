/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include "GraphGenerator.h"

namespace EnsembleClustering {

GraphGenerator::GraphGenerator() {
	// TODO Auto-generated constructor stub

}

GraphGenerator::~GraphGenerator() {
	// TODO Auto-generated destructor stub
}


// TODO: parallel? is insertEdge thread safe?


Graph GraphGenerator::makeErdosRenyiGraph(int64_t n, double p) {
	RandomProbability randP;
	Graph G(n);
	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			if (randP.generate() <= p) {
				G.insertEdge(u, v);
			}
		}
	}
	return G;
}

Graph GraphGenerator::makeRandomGraph(int64_t n, double p) {
	return this->makeErdosRenyiGraph(n, p);	// alias
}

Graph GraphGenerator::makeCircularGraph(int64_t n) {
	// TODO: modernize
	Graph G(n);
	for (int i = 0; i < n; ++i) {
		G.insertEdge(i + 1, ((i+1) % n) + 1);
	}
	return G;
}

Graph GraphGenerator::makeCompleteGraph(int64_t n) {
	// TODO: modernize
	RandomProbability randP;
	Graph G(n);
	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			G.insertEdge(u, v);
		}
	}
	return G;
}



Graph GraphGenerator::makeClusteredRandomGraph(int64_t n, int64_t k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	RandomProbability randP;
	RandomInteger randInt(1, k);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forallNodes([&](node v){
		cluster c = randInt.generate();
		zeta.addToCluster(c, v);
	});

	assert (zeta.numberOfClusters() == k);

	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
				if (randP.generate() <= pin) {
					G.insertEdge(u, v);
				}
			} else {
				if (randP.generate() <= pout) {
					G.insertEdge(u, v);
				}
			}
		}
	}

	return G;
}

} /* namespace EnsembleClustering */
