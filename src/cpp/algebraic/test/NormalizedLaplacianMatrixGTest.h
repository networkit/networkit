/*
 * NormalizedLaplacianMatrixGTest.h
 *
 *  Created on: 25.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NORMALIZEDLAPLACIANMATRIXGTEST_H_
#define NORMALIZEDLAPLACIANMATRIXGTEST_H_

#include "gtest/gtest.h"
#include "../NormalizedLaplacianMatrix.h"
#include "../../graph/Graph.h"

class NormalizedLaplacianMatrixGTest : public testing::Test {
public:
	NormalizedLaplacianMatrixGTest();
	virtual ~NormalizedLaplacianMatrixGTest();
};

#endif /* NORMALIZEDLAPLACIANMATRIXGTEST_H_ */
