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

// the following methdods are commented in order to create linker-errors should they be used before
// they are actually defined (not returning from a function with a returntype != void is UB):

//std::pair<node, node> Sampling::randomEdge(const Graph& G) {
//	// TODO: implement
//}

//node Sampling::randomNeighbor(const Graph& G, node u) {
//	// TODO: implement
//}



} /* namespace NetworKit */

