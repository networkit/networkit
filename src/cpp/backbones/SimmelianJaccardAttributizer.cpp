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

EdgeAttribute SimmelianJaccardAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	EdgeAttribute jaccardAttribute;

	graph.forEdges([&](node u, node v) {
		count maxNeighborhoodSize = std::max(neighbors[u].size(), neighbors[v].size());
		Redundancy redundancy = getOverlap(u, v, neighbors, maxNeighborhoodSize);

		uEdge key = uEdge(u,v);
		if (jaccardAttribute.contains(key))
			jaccardAttribute.set(key, std::min(redundancy.jaccard, jaccardAttribute[key]));
		else
			jaccardAttribute.set(key, redundancy.jaccard);
	});

	return jaccardAttribute;
}

} /* namespace NetworKit */
