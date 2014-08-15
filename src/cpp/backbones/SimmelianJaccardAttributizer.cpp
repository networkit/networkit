/*
 * SimmelianJaccardAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianJaccardAttributizer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

SimmelianJaccardAttributizer::SimmelianJaccardAttributizer() {}

std::vector<double> SimmelianJaccardAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	std::vector<double> jaccardAttribute(graph.upperEdgeIdBound(), 1.0);

	graph.forEdges([&](node u, node v, edgeid eid) {
		count maxNeighborhoodSize = std::max(neighbors[u].size(), neighbors[v].size());
		Redundancy redundancy = getOverlap(u, v, neighbors, maxNeighborhoodSize);

		jaccardAttribute[eid] = std::min(redundancy.jaccard, jaccardAttribute[eid]);
	});

	return jaccardAttribute;
}

} /* namespace NetworKit */
