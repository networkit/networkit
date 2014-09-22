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

Graph GlobalThresholdFilter::calculate(const Graph& graph, const std::vector<double>& attribute) {
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

Graph GlobalThresholdFilter::cloneNodes(const Graph& graph, bool weighted) {
	Graph backboneGraph (graph.upperNodeIdBound(), weighted, false);

	for (node i = 0; i < graph.upperNodeIdBound(); i++) {
		if (!graph.hasNode(i))
			backboneGraph.removeNode(i);
	}
	return backboneGraph;
}

} /* namespace NetworKit */
