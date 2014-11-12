/*
 *
 */

#include <cmath>
#include "GeometricAverageAttributizer.h"
#include "../auxiliary/Log.h"

std::vector< double > NetworKit::GeometricAverageAttributizer::getAttribute(const NetworKit::Graph &g, const std::vector< double > &attribute) {
	std::vector<double> result(g.upperEdgeIdBound());
	
	std::vector<double> nodeSum(g.upperNodeIdBound());
	
	g.parallelForNodes([&](node u) {
		g.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			nodeSum[u] += attribute[eid];
		});
	});
	
	g.parallelForEdges([&](node u, node v, edgeid eid) {
		if (attribute[eid] > 0) {
			result[eid] = attribute[eid] * 1.0 / std::sqrt(nodeSum[u] * nodeSum[v]);
			if (std::isnan(result[eid])) {
				ERROR("Attribute ", attribute[eid], " couldn't be normalized with sum ", nodeSum[u], " and sum ", nodeSum[v]);
			}
		}
	});
	
	return result;
}
