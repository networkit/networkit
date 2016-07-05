/*
 * ChungLu.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include <numeric>

#include "ChungLuGenerator.h"
#include "../graph/GraphBuilder.h"

namespace NetworKit {

ChungLuGenerator::ChungLuGenerator(const std::vector< NetworKit::count > &degreeSequence) :
		StaticDegreeSequenceGenerator(degreeSequence) {
	sum_deg = std::accumulate(seq.begin(), seq.end(), 0);
	n = (count) seq.size();
}

Graph ChungLuGenerator::generate() {
	GraphBuilder gB(n);

	gB.parallelForNodePairs([&](node u, node v) {
		/* Random number in [0, 1] */
		double randVal = Aux::Random::probability();
		/* Probability of edge (u, v): d(u)*d(v)/sum_deg */
		if (randVal < double(seq[u] * seq[v]) / sum_deg) {
			gB.addHalfOutEdge(u, v);
		}
	});
	return gB.toGraph(true,true);
}

} /* namespace NetworKit */
