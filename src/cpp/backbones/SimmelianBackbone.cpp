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

Graph SimmelianBackbone::calculate(const Graph& g, const count& maxRank, const count& minOverlap) {
	ChibaNishizekiTriangleCounter counter;

	edgeCountMap triangles = counter.triangleCounts(g);
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(g, triangles);

	//Create the backbone graph. Edges will be inserted on the fly.
	Graph backboneGraph (g.upperNodeIdBound());

	//TODO: This seems stupid .. Implement a clone method in Graph instead?
	for (node i = 0; i < g.upperNodeIdBound(); i++) {
		if (!g.hasNode(i))
			backboneGraph.removeNode(i);
	}

	g.forEdges([&](node u, node v) {
		count overlap = getOverlap(neighbors[u], neighbors[v], maxRank);
		if (overlap >= minOverlap)
			backboneGraph.addEdge(u, v);
	});

	return backboneGraph;
}

std::vector<RankedNeighbors> SimmelianBackbone::getRankedNeighborhood(const Graph& g, edgeCountMap& triangles) {
	std::vector<RankedNeighbors> neighbors;
	neighbors.resize(g.upperNodeIdBound());

	g.forNodes([&](node u) {
		//Rank ego's alters from strongly to weakly tied.
		g.forNeighborsOf(u, [&](node v) {
			neighbors[u].push_back(RankedEdge(u, v, triangles[uEdge(u, v)]));
		});
		std::sort(neighbors[u].begin(), neighbors[u].end());

		//Calculate the ranks. TODO: These first two steps could be combined for performance gain.
		count currentRank = 0;	//Rank 1 is considered the best.
		count currentSimmelianness = std::numeric_limits<count>::max();
		for (auto& edge : neighbors[u]) {
			if (edge.simmelianness < currentSimmelianness) {
				currentRank++;
				currentSimmelianness = edge.simmelianness;
			}
			edge.rank = currentRank;
		}
	});

	return neighbors;

}

count SimmelianBackbone::getOverlap(const RankedNeighbors& egoNeighbors, const RankedNeighbors& alterNeighbors, const count& maxRank) {
	count overlap = 0;

	std::vector<RankedEdge>::const_iterator egoIt = egoNeighbors.begin();
	std::vector<RankedEdge>::const_iterator alterIt = alterNeighbors.begin();

	std::set<node> egoNeighborsUnmatched;
	std::set<node> alterNeighborsUnmatched;

	//TODO: parameters...
	/*bool allRanks = false;
	int maxRank = allRanks ? std::max(egoNeighbors.size(), alterNeighbors.size()) : 10;*/

	//TODO: identified (nodes that are incident to an edge)
	//TODO: jaccard index
	for (count rank = 1; rank <= maxRank; rank++) {
		matchNeighbors(egoIt, egoNeighbors, egoNeighborsUnmatched, alterNeighborsUnmatched, rank, overlap);
		matchNeighbors(alterIt, alterNeighbors, alterNeighborsUnmatched, egoNeighborsUnmatched, rank, overlap);
	}

	return overlap;
}

/**
 * Helper function used in getOverlap.
 */
void SimmelianBackbone::matchNeighbors(
	std::vector<RankedEdge>::const_iterator& egoIt,
	const RankedNeighbors& egoNeighbors,
	std::set<node>& egoNeighborsUnmatched,
	std::set<node>& alterNeighborsUnmatched,
	const count& rank,
	count& overlap) {

	for (; egoIt != egoNeighbors.end() && egoIt->rank == rank; egoIt++) {
		node other = egoIt->alter;
		if (alterNeighborsUnmatched.erase(other))
			overlap++;
		else
			egoNeighborsUnmatched.insert(other);
	}
}

} /* namespace NetworKit */
