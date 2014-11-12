/*
 *
 */

#include "LinearizeAttribute.h"
#include <algorithm>
#include <numeric>
#include <tuple>
#include "../auxiliary/Random.h"

NetworKit::LinearizeAttribute::LinearizeAttribute(bool inverse) : inverse(inverse) {

}


std::vector< double > NetworKit::LinearizeAttribute::getAttribute(const NetworKit::Graph &g, const std::vector< double > &attribute) {
	std::vector<double> result(g.upperEdgeIdBound());

	// Special case for m = 1
	if (g.numberOfEdges() == 1) {
		g.forEdges([&](node u, node v, edgeid eid) {
			result[eid] = 0.5;
		});
		return result;
	}

	typedef std::tuple<edgeweight, index, edgeid> edgeTuple;
	std::vector<edgeTuple> sorted(g.upperEdgeIdBound(), std::make_tuple(std::numeric_limits<edgeweight>::max(), std::numeric_limits<index>::max(), none));

	g.parallelForEdges([&](node u, node v, edgeid eid) {
		sorted[eid] = std::make_tuple(attribute[eid], Aux::Random::integer(), eid);
	});

	if (inverse) {
		std::sort(sorted.begin(), sorted.end(), std::greater<edgeTuple>());
	} else {
		std::sort(sorted.begin(), sorted.end());
	}

	#pragma omp parallel for
	for (index pos = 0; pos < g.upperEdgeIdBound(); ++pos) {
		edgeid eid = std::get<2>(sorted[pos]);
		if (eid != none) {
			result[eid] = pos * 1.0 / (g.numberOfEdges()-1); // TODO maybe writing back to sorted and then sorting again by eid could be faster (cache efficiency)?
		}
	}

	return result;
}
