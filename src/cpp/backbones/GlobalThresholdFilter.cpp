/*
 * GlobalThresholdFilter.cpp
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#include "GlobalThresholdFilter.h"

namespace NetworKit {

GlobalThresholdFilter::GlobalThresholdFilter(double threshold, bool above) :
		threshold(threshold), above(above) {}

Graph GlobalThresholdFilter::calculate(const Graph& graph, const EdgeAttribute& attribute) {
	//Create an edge-less backbone graph.
	Graph backboneGraph = cloneNodes(graph, false);

	//Re-add the backbone edges.
	graph.forEdges([&](node u, node v) {
		if ((above && attribute[graph.edgeid(u, v)] >= threshold)
				|| (!above && attribute[graph.edgeid(u, v)] <= threshold)) {
			backboneGraph.addEdge(u, v);
		}
	});

	return backboneGraph;
}

} /* namespace NetworKit */
