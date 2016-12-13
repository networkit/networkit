/*
 * MaxentStressGTest.h
 *
 *  Created on: Apr 19, 2016
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_VIZ_TEST_MAXENTSTRESSGTEST_H_
#define NETWORKIT_CPP_VIZ_TEST_MAXENTSTRESSGTEST_H_

#include "gtest/gtest.h"

#include "../../graph/Graph.h"
#include "../Point.h"
#include <vector>
#include <string>


namespace NetworKit {

class MaxentStressGTest : public testing::Test {
public:
	MaxentStressGTest() = default;
	virtual ~MaxentStressGTest() = default;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_VIZ_TEST_MAXENTSTRESSGTEST_H_ */
