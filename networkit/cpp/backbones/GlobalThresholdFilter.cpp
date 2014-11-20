/*
 * GlobalThresholdFilter.cpp
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#include "GlobalThresholdFilter.h"
#include "../graph/GraphBuilder.h"

namespace NetworKit {

GlobalThresholdFilter::GlobalThresholdFilter(const Graph& graph, const std::vector<double>& attribute, const double threshold, bool above) :
		graph(graph), attribute(attribute), threshold(threshold), above(above) {}

Graph GlobalThresholdFilter::calculate() {
	//Create an edge-less backbone graph.
	GraphBuilder builder(graph.upperNodeIdBound(), false, false, true);

	//Re-add the backbone edges.
	graph.balancedParallelForNodes([&](node u) {
		// add each edge in both directions
		graph.forEdgesOf(u, [&](node u, node v, edgeid eid) {
			if ((above && attribute[eid] >= threshold)
			|| (!above && attribute[eid] <= threshold)) {
				builder.addEdge(u, v);
			}
		});
	});

	Graph backboneGraph = builder.toGraph();
	backboneGraph.parallelForNodes([&](node u) {
		if (!graph.hasNode(u)) {
			backboneGraph.removeNode(u);
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
