/*
 * LocalDegreeScore.cpp
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#include "LocalDegreeScore.h"
#include "LocalSimilarityScore.h"

namespace NetworKit {

LocalDegreeScore::LocalDegreeScore(const Graph& G) : EdgeScore<double>(G) {
}

void LocalDegreeScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	std::vector<double> exponents (G.upperEdgeIdBound(), 0.0);

	G.balancedParallelForNodes([&](node i) {
		count d = G.degree(i);

		/**
		 *  The top d^e edges (sorted by degree)
		 * are to be kept in the graph */

		std::vector<AttributizedEdge<count>> neighbors;
		neighbors.reserve(G.degree(i));
		G.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			neighbors.push_back(AttributizedEdge<count>(i, j, eid, G.degree(j)));
		});
		std::sort(neighbors.begin(), neighbors.end());

		/**
		 * By convention, we want to the edges with highest "similarity" or "cohesion" to have values close to 1,
		 * so we invert the range.
		 */

		count rank = 0;
		count numSame = 1;
		count oldValue = 0; // none of the neighbors will have degree 0, so 0 is a safe start value

		#pragma omp critical
		for (auto neighborEdge : neighbors) {
			if (neighborEdge.value != oldValue) {
				rank += numSame;
				numSame = 1;
			} else {
				++numSame;
			}

			edgeid eid = neighborEdge.eid;

			double e = 1.0; // If the node has only one neighbor, the edge should be kept anyway.
			if (d > 1)
				e = 1.0 - (log(rank) / log(d));

			exponents[eid] = std::max(e, exponents[eid]);
		}

	});

	scoreData = std::move(exponents);
	hasRun = true;
}

double LocalDegreeScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double LocalDegreeScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
