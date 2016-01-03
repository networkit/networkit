/*
 * LocalSimilarityScore.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityScore.h"
#include <math.h> //log
#include <set>

namespace NetworKit {

LocalSimilarityScore::LocalSimilarityScore(const Graph& G, const std::vector<count>& triangles) :
	EdgeScore<double>(G), triangles(triangles) {}

void LocalSimilarityScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	/*
	 * For each edge, we calculate the minimum required sparsification exponent e
	 * such that the edge is contained in the sparse graph.
	 */

	std::vector<double> sparsificationExp(G.upperEdgeIdBound(), 0.0);

	G.balancedParallelForNodes([&](node i) {
		count d = G.degree(i);

		/* The top d^e edges (sorted by similarity)
		 * are to be kept in the graph. */

		std::vector<AttributizedEdge<double>> neighbors;
		neighbors.reserve(G.degree(i));
		G.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			double sim = triangles[eid] * 1.0 / (G.degree(i) + G.degree(j) - triangles[eid]);
			neighbors.emplace_back(i, j, eid, sim);
		});
		std::sort(neighbors.begin(), neighbors.end());

		count rank = 1;

		/**
		 * By convention, we want to the edges with highest "similarity" or "cohesion" to have values close to 1,
		 * so we invert the range.
		 */

		#pragma omp critical
		for(std::vector<AttributizedEdge<double>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
			edgeid eid = it->eid;

			double e = 1.0; //If the node has only one neighbor, the edge will be kept anyway.
			if (d > 1)
				e = 1 - (log(rank) / log(d));

			sparsificationExp[eid] = std::max(e, sparsificationExp[eid]);
			rank++;
		}

	});

	scoreData = std::move(sparsificationExp);
	hasRun = true;
}

double LocalSimilarityScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double LocalSimilarityScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
