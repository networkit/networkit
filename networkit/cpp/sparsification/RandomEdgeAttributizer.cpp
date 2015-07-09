/*
 * RandomEdgeAttributizer.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#include "RandomEdgeAttributizer.h"

namespace NetworKit {

RandomEdgeAttributizer::RandomEdgeAttributizer(const Graph& graph) : graph(graph) {
}

std::vector<double> RandomEdgeAttributizer::getAttribute() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	std::vector<double> randomAttribute(graph.upperEdgeIdBound(), 0.0);

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		//double r = Aux::Random::probability();
		//randomAttribute[eid] = randomness * r + (1 - randomness) * attribute[eid];
		randomAttribute[eid] = Aux::Random::probability();
	});

	return randomAttribute;
}

} /* namespace NetworKit */
