/*
 * SimmelianJaccardAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianJaccardAttributizer.h"

namespace NetworKit {

SimmelianJaccardAttributizer::SimmelianJaccardAttributizer() {}

edgeAttribute SimmelianJaccardAttributizer::getAttribute(const Graph& graph, const edgeAttribute& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	edgeAttribute jaccardAttribute;

	graph.forEdges([&](node u, node v) {
		count maxNeighborhoodSize = std::max(neighbors[u].size(), neighbors[v].size());
		Redundancy redundancy = getOverlap(u, v, neighbors, maxNeighborhoodSize);

		uEdge key = uEdge(u,v);
		if (jaccardAttribute.find(key) != jaccardAttribute.end()) {
			double currentJaccard = jaccardAttribute.at(key);
			jaccardAttribute[key] = std::min(redundancy.jaccard, currentJaccard);
		}
		else
			jaccardAttribute[key] = redundancy.jaccard;
	});

	return jaccardAttribute;
}

} /* namespace NetworKit */
