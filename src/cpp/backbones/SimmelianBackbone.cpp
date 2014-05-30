/*
 * SimmelianBackbone.cpp
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianBackbone.h"
#include "TriangleCounter.h"
#include "ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

Graph SimmelianBackbone::calculate(const Graph& g) {
	ChibaNishizekiTriangleCounter counter;

	edgeCountMap triangles = counter.triangleCounts(g);
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(g, triangles);

	return g;
}

count SimmelianBackbone::test() {
	return 42;
}

std::vector<RankedNeighbors> SimmelianBackbone::getRankedNeighborhood(const Graph& g, edgeCountMap& triangles) {
	std::vector<RankedNeighbors> neighbors;

	g.forNodes([&](node u) {
		//Rank ego's alters from strongly to weakly tied.
		g.forNeighborsOf(u, [&](node v) {
			neighbors[u].push_back(SimmelianTie(uEdge(u, v), triangles[uEdge(u, v)]));
		});
		std::sort(neighbors[u].begin(), neighbors[u].end());
	});

	return neighbors;

}

} /* namespace NetworKit */
