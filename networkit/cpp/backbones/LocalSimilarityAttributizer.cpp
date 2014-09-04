/*
 * LocalSimilarityAttributizer.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityAttributizer.h"
#include "ChibaNishizekiTriangleCounter.h"
#include <math.h> //log
#include <set>

namespace NetworKit {

LocalSimilarityAttributizer::LocalSimilarityAttributizer() {}

std::vector<double> LocalSimilarityAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	//Calculate local similarities (using triangle counts)
	std::vector<double> similarity = getLocalSimilarity(graph);

	/*
	* For each edge, we calculate the minimum required sparsification exponent e
	* such that the edge is contained in the backbone.
	*/
	std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), 1.0);

	graph.forNodes([&](node i) {
		count d = graph.degree(i);

		/* The top d^e edges (sorted by similarity)
		 * are to be kept in the backbone */

		std::vector<AttributizedEdge<double>> neighbors;
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			double sim = similarity[eid];
			neighbors.push_back(AttributizedEdge<double>(i, j, eid, sim));
		});
		std::sort(neighbors.begin(), neighbors.end());

		count rank = 1;
		//Top d^e edges are to be retained in the backbone graph.
		//So we calculate the minimum exponent e for each edge that will keep it in the backbone.
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

std::vector<double> LocalSimilarityAttributizer::getLocalSimilarity(const Graph& graph) {
	ChibaNishizekiTriangleCounter triangleAttributizer;
	std::vector<int> triangles = triangleAttributizer.getAttribute(graph, std::vector<int>(graph.upperEdgeIdBound()));
	std::vector<double> similarity(graph.upperEdgeIdBound(), 0.0);
	graph.forEdges([&](node u, node v, edgeid eid) {
		similarity[eid] = ((double) triangles[eid]) / (graph.degree(u) + graph.degree(v) - triangles[eid]);
	});
	return similarity;
}

} /* namespace NetworKit */
