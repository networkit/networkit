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
	// TODO: assign nodes to affinity blocks of d +1 nodes with degree d
}

void BTERGenerator::phaseOne() {
	// TODO: requires: mapping of nodes to affinity blocks - do intelligent implementation
}

void BTERGenerator::phaseTwo() {
	// TODO: requires: expected excess degree of
}

} /* namespace NetworKit */
