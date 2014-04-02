/*
 * LaplacianMatrixGTest.h
 *
 *  Created on: 25.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LAPLACIANMATRIXGTEST_H_
#define LAPLACIANMATRIXGTEST_H_

#include "gtest/gtest.h"
#include "../LaplacianMatrix.h"
#include "../../graph/Graph.h"

class LaplacianMatrixGTest : public testing::Test {
public:
	LaplacianMatrixGTest();
	virtual ~LaplacianMatrixGTest();
};

#endif /* LAPLACIANMATRIXGTEST_H_ */
