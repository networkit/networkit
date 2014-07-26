/*
 * LocalSimilarityAttributizer.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityAttributizer.h"
#include "../structures/Partition.h"
#include "../community/JaccardMeasure.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

LocalSimilarityAttributizer::LocalSimilarityAttributizer() {}

EdgeAttribute LocalSimilarityAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	//EdgeAttribute temp;

	return attribute;
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
