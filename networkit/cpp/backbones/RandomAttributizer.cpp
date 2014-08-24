/*
 * RondomAttributizer.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#include "RandomAttributizer.h"

namespace NetworKit {

RandomAttributizer::RandomAttributizer() {}

std::vector<double> RandomAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	std::vector<double> randomAttribute(graph.upperEdgeIdBound(), 0.0);

	graph.forEdges([&](node u, node v, edgeid eid) {
		randomAttribute[eid] = Aux::Random::probability();
	});

	return randomAttribute;
}

} /* namespace NetworKit */
