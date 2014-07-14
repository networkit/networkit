/*
 * AuxGTest.h
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef AUXGTEST_H_
#define AUXGTEST_H_

// this define is an obscure fix for std::this_thread::sleep_for to work - the issue is described here: http://stackoverflow.com/questions/4438084/stdthis-threadsleep-for-and-gcc
#define _GLIBCXX_USE_NANOSLEEP 1

#include <gtest/gtest.h>

class AuxGTest: public testing::Test {

};

#endif /* AUXGTEST_H_ */

#endif /*NOGTEST */
