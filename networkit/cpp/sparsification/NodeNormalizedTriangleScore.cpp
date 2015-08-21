/*
 * NodeNormalizedTriangleScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include "NodeNormalizedTriangleScore.h"

namespace NetworKit {

NodeNormalizedTriangleScore::NodeNormalizedTriangleScore(const Graph& G, const std::vector<count>& triangles) : EdgeScore<double>(G), triangles(triangles) {
};

void NodeNormalizedTriangleScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	std::vector<double> averageTrianglesPerNode(G.upperNodeIdBound());

	G.balancedParallelForNodes([&](node u) {
		G.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			averageTrianglesPerNode[u] += triangles[eid];
		});
		
		averageTrianglesPerNode[u] /= G.degree(u);
	});

	scoreData.resize(G.upperEdgeIdBound());

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		if (triangles[eid] == 0) {
			if (std::min(averageTrianglesPerNode[u], averageTrianglesPerNode[v]) == 0) {
				scoreData[eid] = G.numberOfNodes() * 1.0 / (G.degree(u) * G.degree(v)); // no triangles incident to any of the two nodes, weight according to degree
			} else { // each node has at least one edge with triangle count > 0 - we can set the weight of this edge to 0
				scoreData[eid] = 0;
			}
		} else {
			scoreData[eid] = triangles[eid] / std::min(averageTrianglesPerNode[u], averageTrianglesPerNode[v]);
		}
	});
	
	hasRun = true;
}

double NodeNormalizedTriangleScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double NodeNormalizedTriangleScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

}
