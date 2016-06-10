/*
 * CommuteTimeDistanceGTest.h
 *
 *  Created on: Jan 17, 2016
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_CENTRALITY_TEST_COMMUTETIMEDISTANCEGTEST_H_
#define NETWORKIT_CPP_CENTRALITY_TEST_COMMUTETIMEDISTANCEGTEST_H_

#include "gtest/gtest.h"
#include "../CommuteTimeDistance.h"

#include <vector>
#include <string>

namespace NetworKit {

using namespace std;

class CommuteTimeDistanceGTest : public testing::Test {
public:
	CommuteTimeDistanceGTest() = default;
	virtual ~CommuteTimeDistanceGTest() = default;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_CENTRALITY_TEST_COMMUTETIMEDISTANCEGTEST_H_ */
