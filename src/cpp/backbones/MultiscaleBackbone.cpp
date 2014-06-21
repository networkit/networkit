/*
 * MultiscaleBackbone.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#include "MultiscaleBackbone.h"

namespace NetworKit {

MultiscaleBackbone::MultiscaleBackbone(const Graph& g, double a) :
		graph(g), alpha(a) {}

Graph MultiscaleBackbone::calculate() {
	Graph backboneGraph = cloneGraphWithoutEdges(graph);

	//The following vector is used for the _local_ normalization of edgeweights.
	//We use a global vector for performance reasons.
	std::vector<edgeweight> normalizedWeights(graph.upperNodeIdBound());

	graph.forNodes([&](node u) {
		int k = graph.degree(u);

		//Normalize edgeweights
		edgeweight maxEdgeweight = 0.0;
		graph.forNeighborsOf(u, [&](node v) {
			maxEdgeweight = std::max(maxEdgeweight, graph.weight(u, v));
		});
		graph.forNeighborsOf(u, [&](node v) {
			normalizedWeights[v] = graph.weight(u, v) / maxEdgeweight;
		});

		//Filter edges by probability calculation
		graph.forNeighborsOf(u, [&](node v) {
			edgeweight p = normalizedWeights[v];
			double a = 1 - (1 - pow(p - 1, k - 1));
			if (a < alpha)
				backboneGraph.setWeight(u, v, graph.weight(u, v));
		});
	});

	return backboneGraph;
}

} /* namespace NetworKit */
