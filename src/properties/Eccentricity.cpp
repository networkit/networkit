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
	static BFS bfs;
	auto dists = bfs.run(G, u);
	auto max_iter = std::max_element(std::begin(dists), std::end(dists));
	return {distance(std::begin(dists), max_iter), *max_iter}; // TODO: unclear
}


} /* namespace NetworKit */

