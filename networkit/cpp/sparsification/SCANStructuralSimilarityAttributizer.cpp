#include "SCANStructuralSimilarityAttributizer.h"

NetworKit::SCANStructuralSimilarityAttributizer::SCANStructuralSimilarityAttributizer(const NetworKit::Graph &graph, const std::vector< NetworKit::count > &triangles) : graph(graph), triangles(triangles) { }

std::vector< double > NetworKit::SCANStructuralSimilarityAttributizer::getAttribute() {
	std::vector<double> attribute(graph.upperEdgeIdBound());

	if (!graph.hasEdgeIds()) throw std::runtime_error("Error, edges must be indexed");

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		attribute[eid] = (triangles[eid] + 1) * 1.0 / std::sqrt((graph.degree(u) + 1)*(graph.degree(v) + 1));
	});

	return attribute;
}

