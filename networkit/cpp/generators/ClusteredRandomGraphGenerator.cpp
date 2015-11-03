/*
 * ClusteredRandomGraphGenerator.cpp
 *
 *  Created on: 28.02.2014
 *      Author: cls
 */

#include "ClusteredRandomGraphGenerator.h"
#include "../structures/Partition.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

ClusteredRandomGraphGenerator::ClusteredRandomGraphGenerator(count n, count k, double pin, double pout) : n(n), k(k), pin(pin), pout(pout) {
}

Graph ClusteredRandomGraphGenerator::generate() {
	assert(pin >= pout);

	Graph G(n);
	// assign nodes evenly to clusters
	Partition zeta(n);
	zeta.setUpperBound(k);
	G.forNodes([&](node v){
		index c = Aux::Random::integer(k-1);
		zeta.addToSubset(c, v);
	});

	assert (zeta.numberOfSubsets() == k);

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

	this->zeta = std::move(zeta);

	G.shrinkToFit();
	return G;

}

Partition ClusteredRandomGraphGenerator::getCommunities() {
	return zeta;
}


} /* namespace NetworKit */


