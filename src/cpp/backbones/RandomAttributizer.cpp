/*
 * RondomAttributizer.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#include "RandomAttributizer.h"

namespace NetworKit {

RandomAttributizer::RandomAttributizer() {}

EdgeAttribute RandomAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	EdgeAttribute randomAttribute(graph.upperEdgeIdBound(), 0.0);

	graph.forEdges([&](node u, node v, edgeid eid) {
		randomAttribute[eid] = Aux::Random::probability();
	});

	return randomAttribute;
}

} /* namespace NetworKit */
