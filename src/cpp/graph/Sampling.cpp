/*
 * Sampling.cpp
 *
 *  Created on: 17.02.2014
 *      Author: cls
 */

#include "Sampling.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

node Sampling::randomNode(const Graph& G) {
	assert (G.numberOfNodes() > 0);
	node v = none;
	do {
		v = Aux::Random::integer(G.upperNodeIdBound());
	} while (!G.hasNode(v));
	return v;
}

std::pair<node, node> Sampling::randomEdge(const Graph& G) {
	// TODO:
}

node Sampling::randomNeighbor(const Graph& G, node u) {
	// TODO:
}


} /* namespace NetworKit */

