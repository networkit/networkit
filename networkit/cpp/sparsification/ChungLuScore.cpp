/*
 * ChungLuScore.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include "ChungLuScore.h"

namespace NetworKit {

ChungLuScore::ChungLuScore(const Graph& G) : EdgeScore<double>(G) {}

void ChungLuScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	scoreData.resize(G.upperEdgeIdBound());
	
	G.parallelForEdges([&](node u, node v, edgeid eid) {
		scoreData[eid] = 1.0 / (G.degree(u) * G.degree(v));
	});
}

} /* namespace NetworKit */
