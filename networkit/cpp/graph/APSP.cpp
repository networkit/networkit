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

APSP::APSP(const Graph& G) : G(G) {}

void APSP::run() {
	std::vector<edgeweight> distanceVector(G.upperNodeIdBound(), 0.0);
	std::vector<std::vector<edgeweight> > distanceMatrix(G.upperNodeIdBound(), distanceVector);
	G.parallelForNodes([&](node u){
		Dijkstra dijk(G, u);
		dijk.run();
		distanceMatrix[u] = dijk.getDistances();
	});


	distances = distanceMatrix;
}

} /* namespace NetworKit */
