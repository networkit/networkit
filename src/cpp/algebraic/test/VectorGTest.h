/*
 * VectorGTest.h
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NOGTEST

#ifndef VECTORGTEST_H_
#define VECTORGTEST_H_

#include <gtest/gtest.h>
#include "../Vector.h"
#include "../../auxiliary/Log.h"
#include <cmath>

class VectorGTest : public testing::Test {
public:
	VectorGTest();
	virtual ~VectorGTest();
};

#endif /* VECTORGTEST_H_ */

#endif
