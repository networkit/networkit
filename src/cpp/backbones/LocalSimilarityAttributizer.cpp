/*
 * LocalSimilarityAttributizer.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityAttributizer.h"
#include "../structures/Partition.h"
#include "../community/JaccardMeasure.h"

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
	//TODO: The following implementation is quite inefficient....
	Graph _graph = graph;

	Partition uNeighbors(1);
	graph.forNeighborsOf(u, [&](node n) {
		uNeighbors.addToSubset(1, n);
	});
	Partition vNeighbors(1);
	graph.forNeighborsOf(v, [&](node n) {
		vNeighbors.addToSubset(1, n);
	});

	JaccardMeasure jaccard;
	return jaccard.getDissimilarity(_graph, uNeighbors, vNeighbors);
}

} /* namespace NetworKit */
