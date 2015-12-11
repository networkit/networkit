/*
 * DenseMatrixGTest.h
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_TEST_DENSEMATRIXGTEST_H_
#define NETWORKIT_CPP_ALGEBRAIC_TEST_DENSEMATRIXGTEST_H_

#include "gtest/gtest.h"

#include "../DenseMatrix.h"
#include "../Vector.h"

namespace NetworKit {

class DenseMatrixGTest : public testing::Test {
	DenseMatrixGTest();
	virtual ~DenseMatrixGTest();
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_TEST_DENSEMATRIXGTEST_H_ */
