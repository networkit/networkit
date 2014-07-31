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
	graph.forEdges([&](node u, node v, edgeid eid) {
		if ((above && attribute[eid] >= threshold)
				|| (!above && attribute[eid] <= threshold)) {
			backboneGraph.addEdge(u, v);
		}
	});

	return backboneGraph;
}

} /* namespace NetworKit */
