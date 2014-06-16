/*
 * GraphGeneratorGTest.h
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef GRAPHGENERATORGTEST_H_
#define GRAPHGENERATORGTEST_H_

#include <gtest/gtest.h>

#include "../GraphGenerator.h"
#include "../../io/GraphIO.h"

namespace NetworKit {

class GraphGeneratorGTest: public testing::Test {
public:
	GraphGeneratorGTest();
	virtual ~GraphGeneratorGTest();
};

} /* namespace NetworKit */
#endif /* GRAPHGENERATORGTEST_H_ */

#endif /*NOGTEST */
