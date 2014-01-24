/*
 * DynCommunityDetector.cpp
 *
 *  Created on: 24.12.2013
 *      Author: cls
 */

#include "DynCommunityDetector.h"

namespace NetworKit {

DynCommunityDetector::DynCommunityDetector() : G(NULL) {
}

void DynCommunityDetector::attachGraph(Graph& G) {
	if (G.numberOfNodes() != 0) {
		throw std::runtime_error("must be initialized with an empty graph");
	}
	this->G = &G;
}


} /* namespace NetworKit */


