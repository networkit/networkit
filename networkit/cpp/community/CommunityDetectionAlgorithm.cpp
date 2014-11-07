/*
 * Clusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "CommunityDetectionAlgorithm.h"

namespace NetworKit {

CommunityDetectionAlgorithm::CommunityDetectionAlgorithm(const Graph& G) : G(G), result(0), hasRun(false) {
}

Partition CommunityDetectionAlgorithm::getPartition() {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return result;
}

std::string CommunityDetectionAlgorithm::toString() const {
	return "TODO: string representation of clusterer";
}

} /* namespace NetworKit */
