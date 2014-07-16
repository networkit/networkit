/*
 * BackboneCalculator.cpp
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#include "BackboneCalculator.h"

namespace NetworKit {

Graph BackboneCalculator::cloneNodes(const Graph& graph, bool weighted) {
	Graph backboneGraph (graph.upperNodeIdBound(), weighted, false);

	for (node i = 0; i < graph.upperNodeIdBound(); i++) {
		if (!graph.hasNode(i))
			backboneGraph.removeNode(i);
	}
	return backboneGraph;
}

} /* namespace NetworKit */
