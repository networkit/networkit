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

	graph.forNodes([&](node u) {
		int k = graph.degree(u);

		graph.forNeighborsOf(u, [&](node v) {
			edgeweight p = graph.weight(u, v);
			double a = 1 - (1 - pow(p - 1, k - 1));
			if (a < alpha)
				backboneGraph.setWeight(u, v, p);
		});
	});

	return backboneGraph;
}

} /* namespace NetworKit */
