/*
 * RandMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "NodeStructuralRandMeasure.h"

namespace NetworKit {

NodeStructuralRandMeasure::~NodeStructuralRandMeasure() {
}

double NodeStructuralRandMeasure::getDissimilarity(Graph& G, Partition& first, Partition& second) {

	count n = G.numberOfNodes();
	assert (n > 0);

	count s11 = 0; 	// number of node pairs for which clusterings aggree
	count s00 = 0;	// number of node pairs for which clusterings disagree

	G.forNodePairs([&](node u, node v){
		if ((first[u] == first[v]) && (second[u] == second[v])) {
			s11 += 1;
		} else if ((first[u] != first[v]) && (second[u] != second[v])) {
			s00 += 1;
		}
	});

	double rand = 1 - ((2 * (s11 + s00)) / (n * (n-1)));

	// assert range [0, 1]
	assert (rand <= 1.0);
	assert (rand >= 0.0);
	return rand;
}

} /* namespace NetworKit */
