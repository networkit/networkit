/*
 *
 */

#include "ChanceCorrectedTriangleAttributizer.h"

std::vector< double > NetworKit::ChanceCorrectedTriangleAttributizer::getAttribute(const NetworKit::Graph &graph, const std::vector< int > &triangles) {
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
