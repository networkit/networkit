/*
 * NodeNormalizedTriangleAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include "NodeNormalizedTriangleAttributizer.h"

namespace NetworKit {

NodeNormalizedTriangleAttributizer::NodeNormalizedTriangleAttributizer(const Graph& graph, const std::vector<int>& triangles) : graph(graph), triangles(triangles) {
};

std::vector< double > NodeNormalizedTriangleAttributizer::getAttribute() {
	std::vector<double> averageTrianglesPerNode(graph.upperNodeIdBound());

	graph.balancedParallelForNodes([&](node u) {
		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			averageTrianglesPerNode[u] += triangles[eid];
		});
		
		averageTrianglesPerNode[u] /= graph.degree(u);
	});

	std::vector<double> result(graph.upperEdgeIdBound());

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		if (triangles[eid] == 0) {
			if (std::min(averageTrianglesPerNode[u], averageTrianglesPerNode[v]) == 0) {
				result[eid] = graph.numberOfNodes() * 1.0 / (graph.degree(u) * graph.degree(v)); // no triangles incident to any of the two nodes, weight according to degree
			} else { // each node has at least one edge with triangle count > 0 - we can set the weight of this edge to 0
				result[eid] = 0;
			}
		} else {
			result[eid] = triangles[eid] / std::min(averageTrianglesPerNode[u], averageTrianglesPerNode[v]);
		}
	});
	
	return result;
}

}
