/*
 * OverlapGTest.h
 *
 *  Created on: 21.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef OVERLAPGTEST_H_
#define OVERLAPGTEST_H_

#include <gtest/gtest.h>
#include <functional>

#include "../RegionGrowingOverlapper.h"
#include "../HashingOverlapper.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/base/ClusteringGenerator.h"


namespace EnsembleClustering {

class OverlapGTest: public testing::Test {

};




} /* namespace EnsembleClustering */
#endif /* OVERLAPGTEST_H_ */
