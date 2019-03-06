/*
 * NeighborhoodDistanceIndex.cpp
 *
 *  Created on: 24.06.2013
 *      Authors: cls, Kolja Esders
 */

#include "../../include/networkit/linkprediction/NeighborhoodDistanceIndex.hpp"
#include "../../include/networkit/linkprediction/NeighborhoodUtility.hpp"

namespace NetworKit {

double NeighborhoodDistanceIndex::runImpl(node u, node v) {
	count uNeighborhood = G->degree(u);
	count vNeighborhood = G->degree(v);
	count intersection = NeighborhoodUtility::getCommonNeighbors(*G, u, v).size();
	return ((double)intersection) / (sqrt(uNeighborhood * vNeighborhood));
}

} /* namespace NetworKit */
