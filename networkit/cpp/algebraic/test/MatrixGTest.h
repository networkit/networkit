/*
 * MatrixGTest.h
 *
 *  Created on: 16.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NOGTEST

#ifndef MATRIXGTEST_H_
#define MATRIXGTEST_H_

#include "gtest/gtest.h"
#include "../Matrix.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../AdjacencyMatrix.h"
#include <math.h>
#include <vector>
#include <utility>


namespace NetworKit {

class MatrixGTest : public testing::Test {
public:
	MatrixGTest();
	virtual ~MatrixGTest();
};


} /* namespace NetworKit */

#endif /* MATRIXGTEST_H_ */

#endif
