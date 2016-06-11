/*
 * ParametrizedGTest.h
 *
 *  Created on: 16.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST


#ifndef PARAMETRIZEDGTEST_H_
#define PARAMETRIZEDGTEST_H_

#include <gtest/gtest.h>


class ParametrizedGTest: public testing::TestWithParam<int> {
};

#endif /* PARAMETRIZEDGTEST_H_ */

#endif /* NOGTEST */
