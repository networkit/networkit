/*
 *
 */

#include "JaccardSimilarityAttributizer.h"

std::vector< double > NetworKit::JaccardSimilarityAttributizer::getAttribute(const NetworKit::Graph &g, const std::vector< int > &triangles) {
	std::vector< double> result(g.upperEdgeIdBound());

	g.parallelForEdges([&](node u, node v, edgeid eid) {
		result[eid] = triangles[eid] * 1.0 / (g.degree(u) + g.degree(v) - triangles[eid]);
	});

	return result;
}
