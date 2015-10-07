/*
 * RandomEdgeScore.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#include "RandomEdgeScore.h"

namespace NetworKit {

RandomEdgeScore::RandomEdgeScore(const Graph& G) : EdgeScore<double>(G) {
}

void RandomEdgeScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}
	scoreData.resize(G.upperEdgeIdBound(), 0.0);

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		//double r = Aux::Random::probability();
		//scoreData[eid] = randomness * r + (1 - randomness) * attribute[eid];
		scoreData[eid] = Aux::Random::probability();
	});
	hasRun = true;
}

double RandomEdgeScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double RandomEdgeScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
