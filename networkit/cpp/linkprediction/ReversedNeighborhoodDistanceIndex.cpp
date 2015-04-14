/*
 * ReversedNeighborhoodDistanceIndex.cpp
 *
 *  Created on: 24.06.2013
 *      Authors: cls, Kolja Esders
 */

#include "ReversedNeighborhoodDistanceIndex.h"
#include "NeighborhoodUtility.h"

namespace NetworKit {

double ReversedNeighborhoodDistanceIndex::runImpl(node u, node v) {
	count uNeighborhood = G->degree(u);
	count vNeighborhood = G->degree(v);
	count intersection = NeighborhoodUtility::getCommonNeighbors(*G, u, v).size();
	if (!G->hasEdge(u, u)) {
		uNeighborhood++;
	}
	if (!G->hasEdge(v, v)) {
		vNeighborhood++;
	}
	return -1 * (1 - ((double) (intersection + 2)) / (sqrt(uNeighborhood * vNeighborhood)));
}

} /* namespace NetworKit */
