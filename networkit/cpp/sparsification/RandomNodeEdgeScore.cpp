/*
 * RandomNodeEdgeScore.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include "RandomNodeEdgeScore.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

RandomNodeEdgeScore::RandomNodeEdgeScore(const Graph& G, double rneRatio) : EdgeScore<double>(G), rneRatio(rneRatio) {
}

void RandomNodeEdgeScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	Graph backbone = G;
	std::vector<double> workScores(G.upperEdgeIdBound(), 0);
	count numRemoved = 0;
	std::vector< std::pair<node, node> > uniformlyRandomEdges;

	while (backbone.numberOfEdges() > 0) {
		if (Aux::Random::real() >= rneRatio) { // uniformly random
			bool edgeFound = false;

			while (!edgeFound) {
				if (uniformlyRandomEdges.empty()) {
					uniformlyRandomEdges = backbone.randomEdges(backbone.numberOfEdges() * (1.0 - rneRatio) + 20);
				}

				auto edge = uniformlyRandomEdges.back();
				uniformlyRandomEdges.pop_back();

				if (backbone.hasEdge(edge.first, edge.second)) {
					edgeid id = backbone.edgeId(edge.first, edge.second);

					workScores[id] = numRemoved * 1.0 / G.numberOfEdges();

					backbone.removeEdge(edge.first, edge.second);

					edgeFound = true;
					++numRemoved;
				}
			}
		} else { // random node - edge
			auto edge = backbone.randomEdge();

			edgeid id = backbone.edgeId(edge.first, edge.second);

			workScores[id] = numRemoved * 1.0 / G.numberOfEdges();

			backbone.removeEdge(edge.first, edge.second);

			++numRemoved;
		}
	}

	scoreData = std::move(workScores);
	hasRun = true;
}

} /* namespace NetworKit */
