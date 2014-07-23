/*
 * SimmelianOverlapAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianOverlapAttributizer.h"
#include <limits>

namespace NetworKit {

SimmelianOverlapAttributizer::SimmelianOverlapAttributizer(count maxRank) :
		maxRank(maxRank) {}

edgeAttribute SimmelianOverlapAttributizer::getAttribute(const Graph& graph, const edgeAttribute& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	edgeAttribute overlapAttribute;

	graph.forEdges([&](node u, node v) {
		Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);

		uEdge key = uEdge(u,v);
		if (overlapAttribute.find(key) != overlapAttribute.end()) {
			double currentOverlap = overlapAttribute.at(key);
			overlapAttribute[key] = std::min((double) redundancy.overlap, currentOverlap);
		}
		else
			overlapAttribute[key] = (double) redundancy.overlap;
	});

	return overlapAttribute;
}

} /* namespace NetworKit */
