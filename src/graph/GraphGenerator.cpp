/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "GraphGenerator.h"

#include "../auxiliary/Random.h"

namespace NetworKit {

GraphGenerator::GraphGenerator() {
	// TODO Auto-generated constructor stub

}

GraphGenerator::~GraphGenerator() {
	// TODO Auto-generated destructor stub
}


// TODO: parallel? is insertEdge thread safe?


Graph GraphGenerator::makeErdosRenyiGraph(count n, double p) {
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		if (Aux::Random::probability() <= p) {
			G.addEdge(u, v);
		}
	});

	return G;
}

Graph GraphGenerator::makeRandomGraph(count n, double p) {
	return this->makeErdosRenyiGraph(n, p);	// alias
}

Graph GraphGenerator::makeCircularGraph(count n) {
	Graph G(n);
	G.forNodes([&](node u){
		G.addEdge(u, (u + 1) % n);
	});
	return G;
}

Graph GraphGenerator::makeCompleteGraph(count n) {
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});
	return G;
}



Graph GraphGenerator::makeClusteredRandomGraph(count n, count k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forNodes([&](node v){
		cluster c = Aux::Random::integer(1, k);
		zeta.addToCluster(c, v);
	});

//	assert (zeta.numberOfClusters() == k);

	G.forNodePairs([&](node u, node v){
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			if (Aux::Random::probability() <= pin) {
				G.addEdge(u, v);
			}
		} else {
			if (Aux::Random::probability() <= pout) {
				G.addEdge(u, v);
			}
		}
	});

	return G;
}

std::pair<Graph, Clustering> GraphGenerator::makeClusteredRandomGraphWithReferenceClustering(
		count n, count k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forNodes([&](node v){
		cluster c = Aux::Random::integer(1, k);
		zeta.addToCluster(c, v);
	});

	G.forNodePairs([&](node u, node v){
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			if (Aux::Random::probability() <= pin) {
				G.addEdge(u, v);
			}
		} else {
			if (Aux::Random::probability() <= pout) {
				G.addEdge(u, v);
			}
		}
	});

	return std::make_pair(G, zeta);
}

Graph GraphGenerator::makeClusteredRandomGraph(Clustering& zeta, double pin,
		double pout) {
	assert (pin >= pout);

	count n = zeta.numberOfNodes();
	Graph G(n);

	G.forNodePairs([&](node u, node v){
		if (zeta.inSameCluster(u, v)) {
			if (Aux::Random::probability() <= pin) {
				G.addEdge(u, v);
			}
		} else {
			if (Aux::Random::probability() <= pout) {
				G.addEdge(u, v);
			}
		}
	});

	return G;
}


} /* namespace NetworKit */
