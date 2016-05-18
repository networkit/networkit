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
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	//Create an edge-less graph.
	GraphBuilder builder(graph.upperNodeIdBound(), false);

	//Re-add the edges of the sparsified graph.
	graph.balancedParallelForNodes([&](node u) {
		// add each edge in both directions
		graph.forEdgesOf(u, [&](node u, node v, edgeid eid) {
			if ((above && attribute[eid] >= threshold)
			|| (!above && attribute[eid] <= threshold)) {
				builder.addHalfEdge(u, v);
			}
		});
	});

	Graph sGraph = builder.toGraph(false);
	// WARNING: removeNode() must not be called in parallel (writes on vector<bool> and does non-atomic decrement of number of nodes)!
	sGraph.forNodes([&](node u) {
		if (!graph.hasNode(u)) {
			sGraph.removeNode(u);
		}
	});

	return sGraph;
}

Graph GlobalThresholdFilter::cloneNodes(const Graph& graph, bool weighted) {
	Graph sparsifiedGraph (graph.upperNodeIdBound(), weighted, false);

	for (node i = 0; i < graph.upperNodeIdBound(); i++) {
		if (!graph.hasNode(i))
			sparsifiedGraph.removeNode(i);
	}
	return sparsifiedGraph;
}

} /* namespace NetworKit */
