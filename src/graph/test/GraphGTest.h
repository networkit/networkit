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
#include <unordered_set>

#include "../../auxiliary/Log.h"
#include "../Graph.h"
#include "../Subgraph.h"
#include "../GraphGenerator.h"
#include "../../graph/GraphGenerator.h"



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
