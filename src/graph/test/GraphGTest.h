/*
 * GraphGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef GRAPHGTEST_H_
#define GRAPHGTEST_H_

#include <gtest/gtest.h>

#include "../../auxiliary/Log.h"
#include "../Graph.h"
#include "../GraphGenerator.h"


namespace NetworKit {

class GraphGTest: public testing::Test {

protected:

	GraphGenerator gen;

public:

	virtual void SetUp();

	virtual void TearDown();

};





} /* namespace NetworKit */
#endif /* GRAPHGTEST_H_ */

#endif /*NOGTEST */
