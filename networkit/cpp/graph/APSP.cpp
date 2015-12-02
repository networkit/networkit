/*
 * APSP.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#include "APSP.h"
#include "../auxiliary/Log.h"
#include "Dijkstra.h"

namespace NetworKit {

APSP::APSP(const Graph& G) : Algorithm(), G(G) {}

void APSP::run() {
	std::vector<edgeweight> distanceVector(G.upperNodeIdBound(), 0.0);
	distances.resize(G.upperNodeIdBound(), distanceVector);
	G.parallelForNodes([&](node u){
		Dijkstra dijk(G, u);
		dijk.run();
		distances[u] = dijk.getDistances();
	});
	hasRun = true;
}

std::string NetworKit::APSP::toString() const {
	return "All-Pairs Shortest Path Algorithm";
}

} /* namespace NetworKit */
