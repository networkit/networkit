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

SimmelianBackbone::SimmelianBackbone(const Graph& g, count maxRank, const count minOverlap) :
		graph(g), parameterized(true), maxRank(maxRank), minOverlap(minOverlap) {}

SimmelianBackbone::SimmelianBackbone(const Graph& g, double treshold) :
		graph(g), parameterized(false), jaccardTreshold(treshold) {}

Graph SimmelianBackbone::calculate() {
	ChibaNishizekiTriangleCounter counter;

	edgeCountMap triangles = counter.triangleCounts(graph);
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, triangles);

	//Create the backbone graph.
	Graph backboneGraph = cloneGraphWithoutEdges(graph);

	graph.forEdges([&](node u, node v) {
		//TODO: This should be refactored.
		count actualMaxRank = parameterized ? maxRank : std::max(neighbors[u].size(), neighbors[v].size());
		Redundancy redundancy = getOverlap(neighbors[u], neighbors[v], actualMaxRank);

		//TODO: extract to a lambda function maybe?
		if ((parameterized && redundancy.overlap >= minOverlap)
				|| (!parameterized && redundancy.jaccard >= jaccardTreshold))
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
		count equals = 1;
		for (auto& edge : neighbors[u]) {
			if (edge.simmelianness < currentSimmelianness) {
				currentRank += equals;
				currentSimmelianness = edge.simmelianness;
				equals = 1;
			} else {
				equals++;
			}
			edge.rank = currentRank;
		}
	});

	return neighbors;

}

Redundancy SimmelianBackbone::getOverlap(const RankedNeighbors& egoNeighbors, const RankedNeighbors& alterNeighbors, const count& maxRank) {
	//Initialization of output values
	count overlap = 0;
	double jaccard = 0.0;

	std::vector<RankedEdge>::const_iterator egoIt = egoNeighbors.begin();
	std::vector<RankedEdge>::const_iterator alterIt = alterNeighbors.begin();

	std::set<node> egoNeighborsUnmatched;
	std::set<node> alterNeighborsUnmatched;

	//TODO: parameters...
	/*bool allRanks = false;
	int maxRank = allRanks ? std::max(egoNeighbors.size(), alterNeighbors.size()) : 10;*/

	//TODO: identified (nodes that are incident to an edge)
	for (count rank = 1; rank <= maxRank; rank++) {
		matchNeighbors(egoIt, egoNeighbors, egoNeighborsUnmatched, alterNeighborsUnmatched, rank, overlap);
		matchNeighbors(alterIt, alterNeighbors, alterNeighborsUnmatched, egoNeighborsUnmatched, rank, overlap);

		double currentJaccard = double(overlap) / double(overlap + egoNeighborsUnmatched.size() + alterNeighborsUnmatched.size());
		jaccard = std::max(currentJaccard, jaccard);
	}

	return Redundancy(overlap, jaccard);
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
