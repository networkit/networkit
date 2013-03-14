/*
 * GraphGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHGTEST_H_
#define GRAPHGTEST_H_

#include <gtest/gtest.h>

#include "../../aux/Log.h"
#include "../Graph.h"
#include "../GraphGenerator.h"


namespace EnsembleClustering {

class GraphGTest: public testing::Test {

protected:

	GraphGenerator gen;

public:

	virtual void SetUp();

	virtual void TearDown();

};





} /* namespace EnsembleClustering */
#endif /* GRAPHGTEST_H_ */
