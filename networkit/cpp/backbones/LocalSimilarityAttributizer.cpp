/*
 * LocalSimilarityAttributizer.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityAttributizer.h"
#include <math.h> //log
#include <set>
#include "../auxiliary/Log.h"

namespace NetworKit {

LocalSimilarityAttributizer::LocalSimilarityAttributizer() {}

std::vector<double> LocalSimilarityAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	/*
	 * For each edge, we calculate the minimum required sparsification exponent e
	 * such that the edge is contained in the backbone.
	 */

	std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), 1.0);

	graph.forNodes([&](node i) {
		count d = graph.degree(i);

		/* The top d^e edges (sorted by similarity)
		 * are to be kept in the backbone */

		std::vector<AttributizedEdge> neighbors;
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			double sim = getSimilarity(graph, i, j);
			neighbors.push_back(AttributizedEdge(i, j, eid, sim));
		});
		std::sort(neighbors.begin(), neighbors.end());

		count rank = 1;
		//Top d^e edges are to be retained in the backbone graph.
		//So we calculate the minimum exponent e for each edge that will keep it in the backbone.
		for(std::vector<AttributizedEdge>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
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

/**
 * Returns the similarity between two nodes.
 */
double LocalSimilarityAttributizer::getSimilarity(const Graph& graph, node u, node v) {
	//Use the jaccard measure as similarity measure.
	/* TODO: The following implementation might be quite inefficient....
	 * can the implementation in community/JaccardMeasure be used?
	 */
	std::set<node> uNeighbors;
	graph.forNeighborsOf(u, [&](node n) {
		uNeighbors.insert(n);
	});

	count inUnion = graph.degree(u);
	count inIntersection = 0;

	graph.forNeighborsOf(v, [&](node n) {
		if (uNeighbors.erase(n))
			inIntersection++;
		else
			inUnion++;
	});

	return (double) inIntersection / (double) inUnion;
}

} /* namespace NetworKit */
