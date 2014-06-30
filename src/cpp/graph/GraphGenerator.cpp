/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "GraphGenerator.h"

#include "../auxiliary/Random.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

// TODO: parallel? is insertEdge thread safe?


Graph GraphGenerator::makeErdosRenyiGraph(count n, double p) {
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		if (Aux::Random::probability() <= p) {
			G.addEdge(u, v);
		}
	});
	G.shrinkToFit();
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
	G.shrinkToFit();
	return G;
}

Graph GraphGenerator::makeCompleteGraph(count n) {
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});
	G.shrinkToFit();
	return G;
}



Graph GraphGenerator::makeClusteredRandomGraph(count n, count k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	// assign nodes evenly to clusters
	Partition zeta(n);
	zeta.setUpperBound(k);
	G.forNodes([&](node v){
		index c = Aux::Random::integer(k-1);
		zeta.addToSubset(c, v);
	});

	if (zeta.numberOfSubsets() != k) {
		WARN("random clustering does not contain k=",k," cluster: ",zeta.numberOfSubsets());
	}


//	assert (zeta.numberOfClusters() == k);

	G.forNodePairs([&](node u, node v){
		if (zeta.subsetOf(u) == zeta.subsetOf(v)) {
			if (Aux::Random::probability() <= pin) {
				G.addEdge(u, v);
			}
		} else {
			if (Aux::Random::probability() <= pout) {
				G.addEdge(u, v);
			}
		}
	});

	G.shrinkToFit();
	return G;
}

std::pair<Graph, Partition> GraphGenerator::makeClusteredRandomGraphWithReferenceClustering(
		count n, count k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	// assign nodes evenly to clusters
	Partition zeta(n);
	zeta.setUpperBound(k);
	G.forNodes([&](node v){
		index c = Aux::Random::integer(k-1);
		zeta.addToSubset(c, v);
	});

	G.forNodePairs([&](node u, node v){
		if (zeta.subsetOf(u) == zeta.subsetOf(v)) {
			if (Aux::Random::probability() <= pin) {
				G.addEdge(u, v);
			}
		} else {
			if (Aux::Random::probability() <= pout) {
				G.addEdge(u, v);
			}
		}
	});

	G.shrinkToFit();
	return std::make_pair(G, zeta);
}

Graph GraphGenerator::makeClusteredRandomGraph(Partition& zeta, double pin,
		double pout) {
	assert (pin >= pout);

	count n = zeta.numberOfElements();
	Graph G(n);

	G.forNodePairs([&](node u, node v){
		if (zeta.inSameSubset(u, v)) {
			if (Aux::Random::probability() <= pin) {
				G.addEdge(u, v);
			}
		} else {
			if (Aux::Random::probability() <= pout) {
				G.addEdge(u, v);
			}
		}
	});

	G.shrinkToFit();
	return G;
}


} /* namespace NetworKit */
