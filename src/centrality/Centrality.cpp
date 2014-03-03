/*
 * Centrality.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "Centrality.h"

namespace NetworKit {


Centrality::Centrality(const Graph& G) : G(G) {
}

double Centrality::score(node v) {
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
	return scoreData;
}



} /* namespace NetworKit */

