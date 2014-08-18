/*
 *
 */

#include "LinearizeAttribute.h"
#include <algorithm>
#include <numeric>

std::vector< double > NetworKit::LinearizeAttribute::getAttribute(const NetworKit::Graph &g, const std::vector< double > &attribute) {
	std::vector<double> result(g.upperEdgeIdBound());

	// Special case for m = 1
	if (g.numberOfEdges() == 1) {
		g.forEdges([&](node u, node v, edgeid eid) {
			result[eid] = 0.5;
		});
		return result;
	}

	std::vector<std::pair<edgeweight, edgeid> > sorted(g.upperEdgeIdBound(), std::make_pair(std::numeric_limits<edgeweight>::max(), none));

	g.parallelForEdges([&](node u, node v, edgeid eid) {
		sorted[eid] = std::make_pair(attribute[eid], eid);
	});

	std::sort(sorted.begin(), sorted.end());

	#pragma omp parallel for
	for (index pos = 0; pos < g.upperEdgeIdBound(); ++pos) {
		if (sorted[pos].second != none) {
			result[sorted[pos].second] = pos * 1.0 / (g.numberOfEdges()-1); // TODO maybe writing back to sorted and then sorting again by eid could be faster (cache efficiency)?
		}
	}

	return result;
}
