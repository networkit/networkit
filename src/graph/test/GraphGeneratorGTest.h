/*
 * GraphGeneratorGTest.h
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHGENERATORGTEST_H_
#define GRAPHGENERATORGTEST_H_

#include <gtest/gtest.h>

#include "../GraphGenerator.h"
#include "../../io/GraphIO.h"

namespace EnsembleClustering {

class GraphGeneratorGTest: public testing::Test {
public:
	GraphGeneratorGTest();
	virtual ~GraphGeneratorGTest();
};

} /* namespace EnsembleClustering */
#endif /* GRAPHGENERATORGTEST_H_ */
