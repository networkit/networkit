/*
 * BellmanFordGTest.h
 *
 *  Created on: Jun 6, 2016
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_TEST_BELLMANFORDGTEST_H_
#define NETWORKIT_CPP_ALGEBRAIC_TEST_BELLMANFORDGTEST_H_

#include "gtest/gtest.h"

#include "../../../graph/Graph.h"

namespace NetworKit {

class BellmanFordGTest : public testing::Test {
public:
	BellmanFordGTest() = default;
	virtual ~BellmanFordGTest() = default;

protected:
	std::vector<double> classicBF(const Graph& graph, node s) const;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_TEST_BELLMANFORDGTEST_H_ */
