/*
 * DynamicCommunityDetector.cpp
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#include "DynamicCommunityDetector.h"

namespace NetworKit {

DynamicCommunityDetector::DynamicCommunityDetector(Graph& G) {
	this->G = &G;
}

DynamicCommunityDetector::DynamicCommunityDetector() {
}

DynamicCommunityDetector::~DynamicCommunityDetector() {
	// TODO Auto-generated destructor stub
}

} /* namespace NetworKit */
