/*
 * DynamicCommunityDetector.cpp
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#include "DynamicCommunityDetector.h"

namespace NetworKit {

DynamicCommunityDetector::DynamicCommunityDetector() : G(NULL) {
}


DynamicCommunityDetector::~DynamicCommunityDetector() {
	// TODO Auto-generated destructor stub
}

std::vector<count> DynamicCommunityDetector::getTimerHistory() {
	return this->timerHistory;
}

} /* namespace NetworKit */
