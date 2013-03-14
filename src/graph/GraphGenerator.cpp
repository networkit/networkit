/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
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


Graph GraphGenerator::makeErdosRenyiGraph(count n, double p) {
	Aux::RandomProbability randP;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		if (randP.generate() <= p) {
			G.insertEdge(u, v);
		}
	});

	return G;
}

Graph GraphGenerator::makeRandomGraph(count n, double p) {
	return this->makeErdosRenyiGraph(n, p);	// alias
}

Graph GraphGenerator::makeCircularGraph(count n) {
	// TODO: modernize
	Graph G(n);
	G.forNodes([&](node u){
		G.insertEdge(u, (u + 1) % n);
	});
	return G;
}

Graph GraphGenerator::makeCompleteGraph(count n) {
	// TODO: modernize
	Aux::RandomProbability randP;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.insertEdge(u, v);
	});
	return G;
}



Graph GraphGenerator::makeClusteredRandomGraph(count n, count k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	Aux::RandomProbability randP;
	Aux::RandomInteger randInt(1, k);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forNodes([&](node v){
		cluster c = randInt.generate();
		zeta.addToCluster(c, v);
	});

//	assert (zeta.numberOfClusters() == k);

	G.forNodePairs([&](node u, node v){
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			if (randP.generate() <= pin) {
				G.insertEdge(u, v);
			}
		} else {
			if (randP.generate() <= pout) {
				G.insertEdge(u, v);
			}
		}
	});

	return G;
}

std::pair<Graph, Clustering> GraphGenerator::makeClusteredRandomGraphWithReferenceClustering(
		count n, count k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	Aux::RandomProbability randP;
	Aux::RandomInteger randInt(1, k);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forNodes([&](node v){
		cluster c = randInt.generate();
		zeta.addToCluster(c, v);
	});

	G.forNodePairs([&](node u, node v){
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			if (randP.generate() <= pin) {
				G.insertEdge(u, v);
			}
		} else {
			if (randP.generate() <= pout) {
				G.insertEdge(u, v);
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

	Aux::RandomProbability randP;
	G.forNodePairs([&](node u, node v){
		if (zeta.inSameCluster(u, v)) {
			if (randP.generate() <= pin) {
				G.insertEdge(u, v);
			}
		} else {
			if (randP.generate() <= pout) {
				G.insertEdge(u, v);
			}
		}
	});

	return G;
}

Graph GraphGenerator::makePreferentialAttachmentGraph(count n, count a) {
	// FIXME: infinite loop

	Graph G(n);

	// all nodes need to have at least degree 1 - create a path
	for (node v = 0; v < (n-1); ++v) {
		G.insertEdge(v, (v + 1));
	}

	count m = G.numberOfEdges(); // number of edges
	count r = 0;

	G.forNodes([&](node u) {
		TRACE("connecting node " << u);
		for (count i = 0; i < a; ++i) { // for all k new edges
			TRACE("2m = " << 2 * m);
			Aux::RandomInteger randInt(0, 2*m);	// TODO: n * k instantiations of RandomInteger are inefficient because random device reads from /dev/random
			r = randInt.generate();
			TRACE("r = " << r);
			for (node v = 0; v < n; ++v) {
				if (r <= G.degree(v)) {
					// select v
					G.insertEdge(u, v);
					TRACE("inserting edge (" << u << "," << v << ")");
					m += 1;
					r = 0;
					break;
				} else {
					TRACE("skipping node " << v);
				}
				r -= G.degree(v);
			}
		}
	});


	return G;
}

} /* namespace EnsembleClustering */
