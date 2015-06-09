/*
 * Centrality.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "Centrality.h"

namespace NetworKit {


Centrality::Centrality(const Graph& G, bool normalized, bool computeEdgeCentrality) : Algorithm(), G(G), normalized(normalized), computeEdgeCentrality(computeEdgeCentrality) {
	if (computeEdgeCentrality && !G.hasEdgeIds()) {
		throw std::runtime_error("For edge centralities to be computed, edges must be indexed first: call G.indexEdges()");
	}
}

double Centrality::score(node v) {
	if (!hasRun) throw std::runtime_error("Call run method first");
	return scoreData.at(v);
}

std::vector<std::pair<node, double> > Centrality::ranking() {
	std::vector<std::pair<node, double> > ranking;
	G.forNodes([&](node v){
		ranking.push_back({v, scoreData[v]});
	});
	std::sort(ranking.begin(), ranking.end(), [](std::pair<node, double> x, std::pair<node, double> y) { return x.second > y.second; });
	return ranking;
}

std::vector<double> Centrality::scores() {
	if (!hasRun) throw std::runtime_error("Call run method first");
	return scoreData;
}

std::vector<double> Centrality::edgeScores() {
	return edgeScoreData;
}

double Centrality::maximum() {
	throw std::runtime_error("Not implemented: Compute the maximum centrality score in the respective centrality subclass.");
}



} /* namespace NetworKit */
