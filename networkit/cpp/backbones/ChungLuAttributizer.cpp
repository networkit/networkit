/*
 *
 */

#include "ChungLuAttributizer.h"

std::vector< double > NetworKit::ChungLuAttributizer::getAttribute(const NetworKit::Graph &g, const std::vector< int > &attribute) {
	std::vector< double > result(g.upperEdgeIdBound());
	
	g.parallelForEdges([&](node u, node v, edgeid eid) {
		result[eid] = 1.0 / (g.degree(u) * g.degree(v));
	});
	
	return result;
}
