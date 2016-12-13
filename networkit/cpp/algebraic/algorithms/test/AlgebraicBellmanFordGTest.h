/*
 * AlgebraicBellmanFordGTest.h
 *
 *  Created on: Jun 6, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_TEST_ALGEBRAICBELLMANFORDGTEST_H_
#define NETWORKIT_CPP_ALGEBRAIC_TEST_ALGEBRAICBELLMANFORDGTEST_H_

#include "gtest/gtest.h"

#include "../../../graph/Graph.h"

namespace NetworKit {

class AlgebraicBellmanFordGTest : public testing::Test {
public:
	AlgebraicBellmanFordGTest() = default;
	virtual ~AlgebraicBellmanFordGTest() = default;

protected:
	std::vector<double> classicBF(const Graph& graph, node s) const;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_TEST_ALGEBRAICBELLMANFORDGTEST_H_ */
