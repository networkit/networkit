/*
 * LocalDegreeAttributizer.cpp
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#include "LocalDegreeAttributizer.h"
#include "LocalSimilarityAttributizer.h"

namespace NetworKit {

LocalDegreeAttributizer::LocalDegreeAttributizer() {}

std::vector<double> LocalDegreeAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), 1.0);

	graph.forNodes([&](node i) {
		count d = graph.degree(i);

		/* The top d^e edges (sorted by degree)
		* are to be kept in the backbone */

		std::vector<AttributizedEdge<count>> neighbors;
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			if (graph.degree(j) > d)
				neighbors.push_back(AttributizedEdge<count>(i, j, eid, graph.degree(j)));
		});
		std::sort(neighbors.begin(), neighbors.end());

		count rank = 1;
		//Top d^e edges are to be retained in the backbone graph.
		//So we calculate the minimum exponent e for each edge that will keep it in the backbone.
		for (auto neighborEdge : neighbors) {
			edgeid eid = neighborEdge.eid;

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
