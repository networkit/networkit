/*
 * EdgeAttributeBlender.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "EdgeAttributeLinearizer.h"
#include <algorithm>
#include <numeric>
#include <tuple>
#include "../auxiliary/Random.h"

namespace NetworKit {

EdgeAttributeLinearizer::EdgeAttributeLinearizer(const Graph& graph, const std::vector<double>& attribute, bool inverse) : graph(graph), attribute(attribute), inverse(inverse) {
}


std::vector< double > EdgeAttributeLinearizer::getAttribute() {
	std::vector<double> result(graph.upperEdgeIdBound());

	// Special case for m = 1
	if (graph.numberOfEdges() == 1) {
		graph.forEdges([&](node u, node v, edgeid eid) {
			result[eid] = 0.5;
		});
		return result;
	}

	typedef std::tuple<edgeweight, index, edgeid> edgeTuple;
	std::vector<edgeTuple> sorted(graph.upperEdgeIdBound(), std::make_tuple(std::numeric_limits<edgeweight>::max(), std::numeric_limits<index>::max(), none));

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		sorted[eid] = std::make_tuple(attribute[eid], Aux::Random::integer(), eid);
	});

	if (inverse) {
		std::sort(sorted.begin(), sorted.end(), std::greater<edgeTuple>());
	} else {
		std::sort(sorted.begin(), sorted.end());
	}

	#pragma omp parallel for
	for (index pos = 0; pos < graph.upperEdgeIdBound(); ++pos) {
		edgeid eid = std::get<2>(sorted[pos]);
		if (eid != none) {
			result[eid] = pos * 1.0 / (graph.numberOfEdges()-1); // TODO maybe writing back to sorted and then sorting again by eid could be faster (cache efficiency)?
		}
	}

	return result;
}

} /* namespace NetworKit */
