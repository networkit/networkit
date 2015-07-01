/*
* StochasticBlockmodel.h
*
*  Created on: 13.08.2014
*      Author: Christian Staudt
*/

#include "StochasticBlockmodel.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

StochasticBlockmodel::StochasticBlockmodel(count n, count nBlocks, const std::vector<index>& membership, const std::vector<std::vector<double> >& affinity)
	: n(n), nBlocks(nBlocks), membership(membership), affinity(affinity) {
	if (affinity.size() != nBlocks) {
		throw std::runtime_error("affinity matrix must be of size nBlocks x nBlocks");
	}
	if (membership.size() != n) {
		throw std::runtime_error("membership list must be of size nNodes");
	}
}


Graph StochasticBlockmodel::generate() {
	Graph G(n);

	G.forNodePairs([&](node u, node v) {
		index a = membership.at(u);
		index b = membership.at(v);
		assert (a < nBlocks);
		assert (b < nBlocks);
		double p = affinity.at(a).at(b);
		double r = Aux::Random::real();
		if (r <= p) {
			G.addEdge(u, v);
		}
	});
	return G;
}

} /* namespace NetworKit */
