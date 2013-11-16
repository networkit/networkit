/*
 * ClusteringCoefficientGTest.h
 *
 *  Created on: Nov 16, 2013
 *      Author: lbarth, dweiss
 */

#ifndef NOGTEST

#ifndef CLUSTERINGCOEFFICIENTGTEST_H_
#define CLUSTERINGCOEFFICIENTGTEST_H_

#include <gtest/gtest.h>
#include "../ClusteringCoefficient.h"

namespace NetworKit {

class ClusteringCoefficientGTest: public testing::Test {
public:
	ClusteringCoefficientGTest();
	virtual ~ClusteringCoefficientGTest();
};

} /* namespace NetworKit */
#endif /* CLUSTERINGCOEFFICIENTGTEST_H_ */

#endif /*NOGTEST */


