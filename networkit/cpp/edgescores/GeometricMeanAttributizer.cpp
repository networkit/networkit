/*
 * GeometricMeanAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include <cmath>
#include "GeometricMeanAttributizer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

GeometricMeanAttributizer::GeometricMeanAttributizer(const Graph& graph, const std::vector<double>& attribute): graph(graph), attribute(attribute) {
}

std::vector< double > GeometricMeanAttributizer::getAttribute() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	std::vector<double> result(graph.upperEdgeIdBound());
	
	std::vector<double> nodeSum(graph.upperNodeIdBound());
	
	graph.parallelForNodes([&](node u) {
		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			nodeSum[u] += attribute[eid];
		});
	});
	
	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		if (attribute[eid] > 0) {
			result[eid] = attribute[eid] * 1.0 / std::sqrt(nodeSum[u] * nodeSum[v]);
			if (std::isnan(result[eid])) {
				ERROR("Attribute ", attribute[eid], " couldn't be normalized with sum ", nodeSum[u], " and sum ", nodeSum[v]);
			}
		}
	});
	
	return result;
}

} /* namespace NetworKit */
