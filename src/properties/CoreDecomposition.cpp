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
	BucketList degrees(0);

	G.forNodes([&](node v) {
		degrees.insert(v, G.degree(v));
	});

	index i = 1;
	Graph G2 = G;
	while (G2.numberOfNodes() > 0) {
		// TODO
	}


	return coreness;
}

} /* namespace NetworKit */
