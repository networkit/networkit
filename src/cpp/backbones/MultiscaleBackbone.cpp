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
		count k = graph.degree(u);

		//Normalize edgeweights
		edgeweight sum = 0.0;
		graph.forNeighborsOf(u, [&](node v) {
			sum += graph.weight(u, v);
		});
		graph.forNeighborsOf(u, [&](node v) {
			normalizedWeights[v] = graph.weight(u, v) / sum;
		});

		//Filter edges by probability calculation
		graph.forNeighborsOf(u, [&](node v) {
			//In case d(u) == 1 and d(v) > 1: consider v only.
			if (k > 1 || graph.degree(v) == 1) {
				edgeweight p = normalizedWeights[v];
				double probability = getProbability(k, p);
				if (probability < alpha) {
					backboneGraph.setWeight(u, v, graph.weight(u, v));
				}
			}
		});
	});

	return backboneGraph;
}

double MultiscaleBackbone::getProbability(count degree, edgeweight normalizedWeight) {
	return 1 - (1 - pow(1 - normalizedWeight, degree - 1));
}

} /* namespace NetworKit */
