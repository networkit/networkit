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

LocalSimilarityAttributizer::LocalSimilarityAttributizer(const Graph& graph, const std::vector<count>& triangles) : graph(graph), triangles(triangles) {}

std::vector<double> LocalSimilarityAttributizer::getAttribute() {
	/*
	 * For each edge, we calculate the minimum required sparsification exponent e
	 * such that the edge is contained in the backbone.
	 */

	std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), 1.0);

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
		//Top d^e edges are to be retained in the backbone graph.
		//So we calculate the minimum exponent e for each edge that will keep it in the backbone.
		#pragma omp critical
		for(std::vector<AttributizedEdge<double>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
			edgeid eid = it->eid;

			double e = 0.0;
			if (d > 1) 			//The node has only one neighbor,, so the edge will be kept anyway.
				e = log(rank) / log(d);

			sparsificationExp[eid] = std::min(e, sparsificationExp[eid]);
			rank++;
		}

	});

	return sparsificationExp;
}

} /* namespace NetworKit */
