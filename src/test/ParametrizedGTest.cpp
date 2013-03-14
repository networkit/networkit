/*
 * ParametrizedGTest.cpp
 *
 *  Created on: 16.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ParametrizedGTest.h"


ParametrizedGTest::ParametrizedGTest() {
	// TODO Auto-generated constructor stub

}

ParametrizedGTest::~ParametrizedGTest() {
	// TODO Auto-generated destructor stub
}


TEST_P(ParametrizedGTest, testParameter) {
	int n = GetParam();
	EXPECT_EQ(n, GetParam());
}


INSTANTIATE_TEST_CASE_P(ParametrizedGTestInstance,
						ParametrizedGTest,
                        ::testing::Values(100, 1000));

