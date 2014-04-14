/*
 * Eccentricity.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "Eccentricity.h"
#include "../graph/BFS.h"

namespace NetworKit {

std::pair<node, count> Eccentricity::getValue(const Graph& G, node u) {
	static BFS bfs(G, u);
	bfs.run();
	auto dists = bfs.getDistances();
	auto max_iter = std::max_element(std::begin(dists), std::end(dists));
	return {std::distance(std::begin(dists), max_iter), *max_iter}; // pair.first is argmax node
}


} /* namespace NetworKit */

