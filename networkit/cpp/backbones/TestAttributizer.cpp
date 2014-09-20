/*
 * TestAttributizer.cpp
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#include "TestAttributizer.h"
#include "LocalSimilarityAttributizer.h"

namespace NetworKit {

TestAttributizer::TestAttributizer(count minDegree, double randomness) : minDegree(minDegree), randomness(randomness) {}

std::vector<double> TestAttributizer::getAttribute(const Graph& graph, const std::vector<int>& triangles) {
	/*std::vector<double> test(graph.upperEdgeIdBound(), 0.0);

	double maxt = (double) *std::max_element(std::begin(triangles), std::end(triangles));

	graph.forEdges([&](node u, node v, edgeid eid) {
		if (graph.degree(u) <= minDegree || graph.degree(v) <= minDegree) {
			test[eid] = 1.0;
		}

		double random = Aux::Random::probability();
		test[eid] = (randomness * random) + (1-randomness) * (triangles[eid] / maxt);
	});

	return test;*/

	std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), 1.0);

	graph.forNodes([&](node i) {
		count d = graph.degree(i);

		/* The top d^e edges (sorted by degree)
		* are to be kept in the backbone */

		std::vector<AttributizedEdge<count>> neighbors;
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			neighbors.push_back(AttributizedEdge<count>(i, j, eid, graph.degree(j)));
		});
		std::sort(neighbors.begin(), neighbors.end());

		count rank = 1;
		//Top d^e edges are to be retained in the backbone graph.
		//So we calculate the minimum exponent e for each edge that will keep it in the backbone.
		for(std::vector<AttributizedEdge<count>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
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
