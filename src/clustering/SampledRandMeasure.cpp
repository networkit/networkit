/*
 * SampledRandMeasure.cpp
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#include "SampledRandMeasure.h"

#include "../auxiliary/Random.h"

namespace NetworKit {

SampledRandMeasure::SampledRandMeasure(count maxSamples) : maxSamples(maxSamples) {

}

SampledRandMeasure::~SampledRandMeasure() {
	// TODO Auto-generated destructor stub
}

double SampledRandMeasure::getDissimilarity(Graph& G, Partition& first, Partition& second) {
	assert (G.numberOfNodes() > 0);
	assert (G.numberOfNodes() > maxSamples);

	count n11 = 0; 	// number of node pairs for which clusterings aggree
	count n00 = 0;	// number of node pairs for which clusterings disagree


	count nSamples = 0;

	index z = G.upperNodeIdBound();

	while (nSamples < maxSamples) {
		node u = Aux::Random::integer(z);
		node v = Aux::Random::integer(z);
		if (u != v) { // nodes should be distinct
			if (G.hasNode(u) && G.hasNode(v)) { // nodes should exist in the graph
				if ((first[u] == first[v]) && (second[u] == second[v])) {
					n11 += 1;
				} else if ((first[u] != first[v]) && (second[u] != second[v])) {
					n00 += 1;
				}
				nSamples += 1;
			}
		}
	}

	DEBUG("n11 = " , n11 , " n00 = " , n00 , " nSamples = " , nSamples);

	double dis = 1 - ((n00 + n11) / (double) nSamples);
	return dis;
}

} /* namespace NetworKit */
