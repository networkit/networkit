/*
 * EdgeScoreLinearizer.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "EdgeScoreLinearizer.h"
#include <numeric>
#include <tuple>
#include "../auxiliary/Random.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

EdgeScoreLinearizer::EdgeScoreLinearizer(const Graph& G, const std::vector<double>& attribute, bool inverse) : EdgeScore<double>(G), attribute(attribute), inverse(inverse) {
}


void EdgeScoreLinearizer::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	scoreData.resize(G.upperEdgeIdBound());

	// Special case for m = 1
	if (G.numberOfEdges() == 1) {
		G.forEdges([&](node u, node v, edgeid eid) {
			scoreData[eid] = 0.5;
		});
	} else {
		typedef std::tuple<edgeweight, index, edgeid> edgeTuple;
		std::vector<edgeTuple> sorted(G.upperEdgeIdBound(), std::make_tuple(std::numeric_limits<edgeweight>::max(), std::numeric_limits<index>::max(), none));

		G.parallelForEdges([&](node u, node v, edgeid eid) {
			sorted[eid] = std::make_tuple(attribute[eid], Aux::Random::integer(), eid);
		});

		if (inverse) {
			Aux::Parallel::sort(sorted.begin(), sorted.end(), std::greater<edgeTuple>());
		} else {
			Aux::Parallel::sort(sorted.begin(), sorted.end());
		}

		#pragma omp parallel for
		for (index pos = 0; pos < G.upperEdgeIdBound(); ++pos) {
			edgeid eid = std::get<2>(sorted[pos]);
			if (eid != none) {
				scoreData[eid] = pos * 1.0 / (G.numberOfEdges()-1); // TODO maybe writing back to sorted and then sorting again by eid could be faster (cache efficiency)?
			}
		}
	}

	hasRun = true;
}

double EdgeScoreLinearizer::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double EdgeScoreLinearizer::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
