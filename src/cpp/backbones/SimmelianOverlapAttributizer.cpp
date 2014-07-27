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

EdgeAttribute SimmelianOverlapAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	EdgeAttribute overlapAttribute;

	graph.forEdges([&](node u, node v) {
		Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);

		uEdge key = uEdge(u,v);
		if (overlapAttribute.contains(key))
			overlapAttribute.set(key, std::min((double) redundancy.overlap, overlapAttribute[key]));
		else
			overlapAttribute.set(key, (double) redundancy.overlap);
	});

	return overlapAttribute;
}

} /* namespace NetworKit */
