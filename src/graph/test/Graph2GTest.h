/*
 * Graph2GTest.h
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef GRAPH2GTEST_H_
#define GRAPH2GTEST_H_

#include <gtest/gtest.h>

#include "../Graph.h"

namespace NetworKit {

class Graph2GTest: public testing::Test {
public:
	Graph2GTest();
	virtual ~Graph2GTest();
};

} /* namespace NetworKit */
#endif /* GRAPH2GTEST_H_ */

#endif /*NOGTEST */
