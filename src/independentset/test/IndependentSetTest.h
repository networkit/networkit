/*
 * IndependentSetTest.h
 *
 *  Created on: 27.02.2013
 *      Author: cls
 */

#ifndef INDEPENDENTSETTEST_H_
#define INDEPENDENTSETTEST_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../independentset/Luby.h"

namespace EnsembleClustering {

class IndependentSetTest: public testing::Test {
public:
	IndependentSetTest();
	virtual ~IndependentSetTest();
};

} /* namespace EnsembleClustering */
#endif /* INDEPENDENTSETTEST_H_ */
