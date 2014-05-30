/*
 * SimmelianBackbone.cpp
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianBackbone.h"
#include "TriangleCounter.h"
#include "ChibaNishizekiTriangleCounter.h"
#include <limits>

namespace NetworKit {

Graph SimmelianBackbone::calculate(const Graph& g) {
	ChibaNishizekiTriangleCounter counter;

	edgeCountMap triangles = counter.triangleCounts(g);
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(g, triangles);

	//Create the backbone graph. Edges will be inserted on the fly.
	Graph backboneGraph (g.upperNodeIdBound());

	g.forEdges([&](node u, node v) {
		uEdge edge = uEdge(u, v);

		int overlap = getOverlap(neighbors[u], neighbors[v]);
	});

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
		std::sort(neighbors[u].begin(), neighbors[u].end(), std::greater<SimmelianTie>());

		//Calculate the ranks. TODO: These first two steps could be combined for performance gain.
		int currentRank = 0;
		int currentSimmelianness = std::numeric_limits<int>::max();
		g.forNeighborsOf(u, [&](node v) {
			if (neighbors[u][v].simmelianness < currentSimmelianness) {
				currentRank++;
				currentSimmelianness = neighbors[u][v].simmelianness;
			}
			neighbors[u][v].rank = currentRank;
		});
	});

	return neighbors;

}

int SimmelianBackbone::getOverlap(const RankedNeighbors& first, const RankedNeighbors& second) {
	int overlap = 0;

	for (std::vector<SimmelianTie>::const_iterator firstIt = first.begin(); firstIt != first.end(); firstIt++) {

	}

	return overlap;
}

} /* namespace NetworKit */
