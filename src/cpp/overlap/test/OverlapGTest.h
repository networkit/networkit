/*
 * OverlapGTest.h
 *
 *  Created on: 21.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef OVERLAPGTEST_H_
#define OVERLAPGTEST_H_

#include <gtest/gtest.h>
#include <functional>

#include "../RegionGrowingOverlapper.h"
#include "../HashingOverlapper.h"
#include "../../graph/GraphGenerator.h"
#include "../../community/ClusteringGenerator.h"
#include "../../auxiliary/Debug.h"


namespace NetworKit {

class OverlapGTest: public testing::Test {

};




} /* namespace NetworKit */
#endif /* OVERLAPGTEST_H_ */

#endif /* NOGTEST */
