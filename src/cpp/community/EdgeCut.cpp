/*
 * EdgeCut.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#include "EdgeCut.h"

namespace NetworKit {

EdgeCut::EdgeCut() {

}

EdgeCut::~EdgeCut() {

}

// TODO: unit test
double EdgeCut::getQuality(const Partition& zeta, const Graph& G) {
	double cutWeight = 0.0;
	G.forEdges([&](node u, node v) {
		if (zeta[u] != zeta[v]) {
			cutWeight += G.weight(u, v);
		}
	});
	return cutWeight;
}

} /* namespace NetworKit */
