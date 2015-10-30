/*
 * EdgeScoreBlender.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "EdgeScoreBlender.h"

namespace NetworKit {

EdgeScoreBlender::EdgeScoreBlender(const Graph &G, const std::vector< double > &attribute0, const std::vector< double > &attribute1, const std::vector< bool > &selection) :
EdgeScore(G), attribute0(attribute0), attribute1(attribute1), selection(selection) {}

void EdgeScoreBlender::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	scoreData.resize(G.upperEdgeIdBound());

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		scoreData[eid] = (selection[eid] ? attribute1[eid] : attribute0[eid]);
	});

	hasRun = true;
}

double EdgeScoreBlender::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double EdgeScoreBlender::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
