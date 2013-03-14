/*
 * RandMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandMeasure.h"

namespace EnsembleClustering {

RandMeasure::RandMeasure() {
	// TODO Auto-generated constructor stub

}

RandMeasure::~RandMeasure() {
	// TODO Auto-generated destructor stub
}

double RandMeasure::getDissimilarity(Graph& G, Clustering& first, Clustering& second) {

	int64_t n = G.numberOfNodes();
	assert (n > 0);

	int64_t s11 = 0; 	// number of node pairs for which clusterings aggree
	int64_t s00 = 0;	// number of node pairs for which clusterings disagree

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

} /* namespace EnsembleClustering */
