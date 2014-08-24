/*
 * SimmelianAttributizer.cpp
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianAttributizer.h"
#include <limits>

namespace NetworKit {

std::vector<RankedNeighbors> SimmelianAttributizer::getRankedNeighborhood(const Graph& g, const std::vector<int>& triangles) {
	std::vector<RankedNeighbors> neighbors;
	neighbors.resize(g.upperNodeIdBound());

	g.forNodes([&](node u) {
		//Sort ego's alters from strongly to weakly tied.
		g.forNeighborsOf(u, [&](node _u, node v, edgeid eid) {
			count triangleCount = round(triangles[eid]);
			neighbors[u].push_back(RankedEdge(u, v, triangleCount));
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

Redundancy SimmelianAttributizer::getOverlap(	const node& ego,
												const node& alter,
												const std::vector<RankedNeighbors>& neighbors,
												const count& maxRank) {
	//Initialization of output values
	Redundancy result = Redundancy(0, 0.0);

	std::vector<RankedEdge>::const_iterator egoIt = neighbors[ego].begin();
	std::vector<RankedEdge>::const_iterator alterIt = neighbors[alter].begin();

	std::set<node> egoNeighborsUnmatched;
	std::set<node> alterNeighborsUnmatched;

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
void SimmelianAttributizer::matchNeighbors(
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
