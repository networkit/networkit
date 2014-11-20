/*
 * ChangeCorrectedTriangleAttributizer.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include "ChanceCorrectedTriangleAttributizer.h"

namespace NetworKit {

ChanceCorrectedTriangleAttributizer::ChanceCorrectedTriangleAttributizer(const Graph& graph, const std::vector<count>& triangles) : graph(graph), triangles(triangles) {
}

std::vector< double > ChanceCorrectedTriangleAttributizer::getAttribute() {
	std::vector<double> result(graph.upperEdgeIdBound(), 0);
	
	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		if (triangles[eid] > 0) {
			result[eid] = triangles[eid] * (graph.numberOfNodes() - 2) * 1.0 / ((graph.degree(u) - 1) * (graph.degree(v) - 1));
		} else if (graph.degree(u) == 1 || graph.degree(v) == 1) {
			result[eid] = 1;
		}
	});
	
	return result;
}

} /* namespace NetworKit */
