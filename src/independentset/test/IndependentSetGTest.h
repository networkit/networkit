/*
 * IndependentSetTest.h
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef INDEPENDENTSETGTEST_H_
#define INDEPENDENTSETGTEST_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../independentset/Luby.h"

namespace EnsembleClustering {

class IndependentSetGTest: public testing::Test {
public:
	IndependentSetGTest();
	virtual ~IndependentSetGTest();
};

} /* namespace EnsembleClustering */
#endif /* INDEPENDENTSETTEST_H_ */
