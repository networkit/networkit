/*
 * SimmelianJaccardAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianJaccardAttributizer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

SimmelianJaccardAttributizer::SimmelianJaccardAttributizer(const Graph& graph, const std::vector<int>& triangles) : SimmelianAttributizer(graph, triangles) {
}

std::vector<double> SimmelianJaccardAttributizer::getAttribute() {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, triangles);
	std::vector<double> jaccardAttribute(graph.upperEdgeIdBound(), 1.0);

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		count maxNeighborhoodSize = std::max(neighbors[u].size(), neighbors[v].size());
		Redundancy redundancy = getOverlap(u, v, neighbors, maxNeighborhoodSize);

		jaccardAttribute[eid] = redundancy.jaccard;
	});

	return jaccardAttribute;
}

} /* namespace NetworKit */
