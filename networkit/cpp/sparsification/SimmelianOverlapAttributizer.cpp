/*
 * SimmelianOverlapAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianOverlapAttributizer.h"
#include <limits>

namespace NetworKit {

SimmelianOverlapAttributizer::SimmelianOverlapAttributizer(const Graph& graph, const std::vector<count>& triangles, count maxRank) :
		SimmelianAttributizer(graph, triangles), maxRank(maxRank) {}

std::vector<double> SimmelianOverlapAttributizer::getAttribute() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, triangles);
	std::vector<double> overlapAttribute(graph.upperEdgeIdBound(), 0.0);

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);

		overlapAttribute[eid] = (double) redundancy.overlap;
	});

	return overlapAttribute;
}

} /* namespace NetworKit */
