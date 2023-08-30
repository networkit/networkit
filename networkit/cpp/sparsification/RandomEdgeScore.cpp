/*
 * RandomEdgeScore.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#include <networkit/sparsification/RandomEdgeScore.hpp>

namespace NetworKit {

RandomEdgeScore::RandomEdgeScore(const Graph &G) : EdgeScore<double>(G) {}

void RandomEdgeScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }
    scoreData.resize(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges(
        [&](node, node, edgeid eid) { scoreData[eid] = Aux::Random::probability(); });
    hasRun = true;
}

double RandomEdgeScore::score(node u, node v) {
    if (hasRun) {
        return scoreData[G->edgeId(u, v)];
    }
    throw std::runtime_error("Call run() prior to using score().");
}

double RandomEdgeScore::score(edgeid eid) {
    if (hasRun) {
        return scoreData[eid];
    }
    throw std::runtime_error("Call run() prior to using score().");
}

} /* namespace NetworKit */
