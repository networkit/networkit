/*
 * BTERGenerator.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "BTERGenerator.h"

namespace NetworKit {

BTERGenerator::BTERGenerator(std::vector<count> degreeDistribution,
		std::vector<count> clusteringCoefficients) {
}

BTERGenerator::~BTERGenerator() {
	// TODO Auto-generated destructor stub
}

Graph BTERGenerator::generate() {
}

void BTERGenerator::preprocessing() {

	std::vector<count> degreeDistribution; // TODO: parameter

	double beta; // TODO:

	// TODO: assign nodes to affinity blocks of d +1 nodes with degree d

	// TODO: precompute index i_d for first node of each degree d
	// TODO: compute number of nodes n'_d with degree greater d

	count dMax; // maximum degree

	// number of nodes from least degree to gretes, except degree-1 nodes are last
	std::vector<index> firstNodeWithDegree;
	firstNodeWithDegree[2] = 1; // TODO: should indices be zero-based?
	for (count d = 3; d <= dMax; d++) {
		firstNodeWithDegree[d] = firstNodeWithDegree[d - 1] + degreeDistribution[d - 1];
	}
	firstNodeWithDegree[1] = firstNodeWithDegree[dMax] + degreeDistribution[dMax];

	// compute number of nodes with degree greater than d
	std::vector<count> higherDegreeNodes;
	for (count d = 1; d < dMax; d++) {
		for (count d_ = d + 1; d_ <= dMax; d_++) {
			higherDegreeNodes[d] += degreeDistribution[d_]; // TODO: check that higherDegreeNodes[dMax]Ê== 0
		}
	}

	// handle degree-1 nodes
	count n1Fill = beta * degreeDistribution[1];
	double w1 = 0.5 * n1Fill;
	count r1Fill = 1;

	// main loop
	count g = 0; // ?
	count nStarFill = 0; // ?
	count dStar = 0; // ?

	for (count d = 2; d <= dMax; d++) {
		if (nStarFill > 0) {

		}
	}

}

void BTERGenerator::phaseOne() {
	// TODO: requires: mapping of nodes to affinity blocks - do intelligent implementation
}

void BTERGenerator::phaseTwo() {
	// TODO: requires: expected excess degree of
}

} /* namespace NetworKit */
