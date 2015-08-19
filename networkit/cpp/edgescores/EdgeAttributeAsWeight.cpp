/*
 * EdgeAttributeAsWeight.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "EdgeAttributeAsWeight.h"

namespace NetworKit {

EdgeAttributeAsWeight::EdgeAttributeAsWeight(const Graph& graph, const std::vector<double>& attribute, bool squared, edgeweight offset, edgeweight factor) :
		graph(graph), attribute(attribute), squared(squared), offset(offset), factor(factor) {
}

Graph EdgeAttributeAsWeight::calculate() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	Graph result(graph, true, false);

	if (squared) {
		graph.parallelForEdges([&](node u, node v, edgeid eid) {
			result.setWeight(u, v, offset + factor * attribute[eid] * attribute[eid]);
		});
	} else {
		graph.parallelForEdges([&](node u, node v, edgeid eid) {
			result.setWeight(u, v, offset + factor * attribute[eid]);
		});
	}

	return result;
}

} // namespace NetworKit
