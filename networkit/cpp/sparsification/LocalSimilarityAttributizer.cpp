/*
 * LocalSimilarityAttributizer.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityAttributizer.h"
#include <math.h> //log
#include <set>

namespace NetworKit {

LocalSimilarityAttributizer::LocalSimilarityAttributizer(const Graph& graph, const std::vector<count>& triangles) :
	graph(graph), triangles(triangles) {}

std::vector<double> LocalSimilarityAttributizer::getAttribute() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	/*
	 * For each edge, we calculate the minimum required sparsification exponent e
	 * such that the edge is contained in the backbone.
	 */

	std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), 0.0);

	graph.balancedParallelForNodes([&](node i) {
		count d = graph.degree(i);

		/* The top d^e edges (sorted by similarity)
		 * are to be kept in the backbone */

		std::vector<AttributizedEdge<double>> neighbors;
		neighbors.reserve(graph.degree(i));
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			double sim = triangles[eid] * 1.0 / (graph.degree(i) + graph.degree(j) - triangles[eid]);
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

	return sparsificationExp;
}

} /* namespace NetworKit */
