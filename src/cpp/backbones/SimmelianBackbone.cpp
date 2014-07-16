/*
 * SimmelianBackbone.cpp
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianBackbone.h"
#include "ChibaNishizekiTriangleCounter.h"
#include <limits>

namespace NetworKit {

SimmelianBackbone::SimmelianBackbone(count maxRank, const count minOverlap) :
		parameterized(true), maxRank(maxRank), minOverlap(minOverlap) {}

SimmelianBackbone::SimmelianBackbone(double treshold) :
		parameterized(false), jaccardTreshold(treshold) {}

Graph SimmelianBackbone::calculate(const Graph& graph) {
	ChibaNishizekiTriangleCounter counter;

	edgeCountMap triangles = counter.triangleCounts(graph);
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, triangles);

	//Create an edge-less backbone graph.
	Graph backboneGraph = cloneNodes(graph, false);

	//Re-add the backbone edges.
	if (parameterized) {
		graph.forEdges([&](node u, node v) {
			Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);
			if (redundancy.overlap >= minOverlap)
				backboneGraph.addEdge(u, v);
		});
	} else {
		graph.forEdges([&](node u, node v) {
			count maxNeighborhoodSize = std::max(neighbors[u].size(), neighbors[v].size());
			Redundancy redundancy = getOverlap(u, v, neighbors, maxNeighborhoodSize);

			if (redundancy.jaccard >= jaccardTreshold)
				backboneGraph.addEdge(u, v);
		});
	}

	return backboneGraph;
}

std::vector<RankedNeighbors> SimmelianBackbone::getRankedNeighborhood(const Graph& g, edgeCountMap& triangles) {
	std::vector<RankedNeighbors> neighbors;
	neighbors.resize(g.upperNodeIdBound());

	g.forNodes([&](node u) {
		//Sort ego's alters from strongly to weakly tied.
		g.forNeighborsOf(u, [&](node v) {
			neighbors[u].push_back(RankedEdge(u, v, triangles[uEdge(u, v)]));
		});
		std::sort(neighbors[u].begin(), neighbors[u].end());

		//Calculate the ranks.
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

Redundancy SimmelianBackbone::getOverlap(	const node& ego,
											const node& alter,
											const std::vector<RankedNeighbors> neighbors,
											const count& maxRank) {
	//Initialization of output values
	Redundancy result = Redundancy(0, 0.0);

	std::vector<RankedEdge>::const_iterator egoIt = neighbors[ego].begin();
	std::vector<RankedEdge>::const_iterator alterIt = neighbors[alter].begin();

	std::set<node> egoNeighborsUnmatched;
	std::set<node> alterNeighborsUnmatched;

	//TODO: identified parameter? (nodes that are incident to an edge)
	for (count rank = 1; rank <= maxRank; rank++) {
		matchNeighbors(ego, alter, true, egoIt, neighbors[ego], egoNeighborsUnmatched, alterNeighborsUnmatched, rank, result.overlap);
		matchNeighbors(alter, ego, false, alterIt, neighbors[alter], alterNeighborsUnmatched, egoNeighborsUnmatched, rank, result.overlap);

		double currentJaccard = double(result.overlap) /
				double(result.overlap + egoNeighborsUnmatched.size() + alterNeighborsUnmatched.size());

		result.jaccard = std::max(currentJaccard, result.jaccard);
	}

	return result;
}

/**
 * Helper function used in getOverlap. Adds the intersection of
 * egoNeighbors and alterNeighborsUnmatched to overlap.
 */
void SimmelianBackbone::matchNeighbors(
	const node& ego,
	const node& alter,
	const bool& reciprocityCheck,
	std::vector<RankedEdge>::const_iterator& egoIt,
	const RankedNeighbors& egoNeighbors,
	std::set<node>& egoNeighborsUnmatched,
	std::set<node>& alterNeighborsUnmatched,
	const count& rank,
	count& overlap) {

	for (; egoIt != egoNeighbors.end() && egoIt->rank == rank; egoIt++) {
		node other = egoIt->alter;

		//We count nodes that are incident to the currently considered edge, as well.
		if (reciprocityCheck && other == alter)
			other = ego;

		if (alterNeighborsUnmatched.erase(other))
			overlap++;
		else
			egoNeighborsUnmatched.insert(other);
	}
}

} /* namespace NetworKit */
