/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition() {
	// TODO Auto-generated constructor stub

}

CoreDecomposition::~CoreDecomposition() {
	// TODO Auto-generated destructor stub
}

std::vector<count> CoreDecomposition::run(const Graph& G) {
	std::vector<count> coreness;

	G.forNodes([&](node v) {
		// TODO: fill data structure
	});

	index i = 1;
	Graph G2 = G;
	while (G2.numberOfNodes() > 0) {
		// TODO: main loop
	}

	return coreness;
}

} /* namespace NetworKit */
