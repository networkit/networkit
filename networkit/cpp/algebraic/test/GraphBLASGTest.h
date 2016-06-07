/*
 * GraphBLASGTest.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_TEST_GRAPHBLASGTEST_H_
#define NETWORKIT_CPP_ALGEBRAIC_TEST_GRAPHBLASGTEST_H_

#include "../GraphBLAS.h"

#include "gtest/gtest.h"

namespace NetworKit {

class GraphBLASGTest : public testing::Test {
public:
	GraphBLASGTest() = default;
	virtual ~GraphBLASGTest() = default;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_TEST_GRAPHBLASGTEST_H_ */
