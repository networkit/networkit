/*
 * NeighborhoodDistanceIndex.cpp
 *
 *  Created on: 24.06.2013
 *      Authors: cls, Kolja Esders
 */

#include "NeighborhoodDistanceIndex.h"

namespace NetworKit {

NeighborhoodDistanceIndex::NeighborhoodDistanceIndex(const Graph& G) : LinkPredictor(G) {
}

double NeighborhoodDistanceIndex::runImpl(node u, node v) {
	count inter = 0;
	int neighborhood1 = G->degree(u);
	int neighborhood2 = G->degree(v);

	G->forNeighborsOf(u, [&](node x){
		if (x != v && x!= u ) {
			if (G->hasEdge(x, v)) {
				inter++;
			}
		}
	});
	if (!G->hasEdge(u, u)) {
		neighborhood1++;
	}
	if (!G->hasEdge(v, v)) {
		neighborhood2++;
	}

	return (1 - ((double) (inter + 2)) / (sqrt(neighborhood2 * neighborhood1)));
}

} /* namespace NetworKit */
