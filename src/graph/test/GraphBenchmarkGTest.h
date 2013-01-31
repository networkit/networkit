/*
 * GraphBenchmarkGTest.h
 *
 *  Created on: 30.01.2013
 *      Author: cls
 */

#ifndef GRAPHBENCHMARKGTEST_H_
#define GRAPHBENCHMARKGTEST_H_

#include <gtest/gtest.h>

#include "../GraphGenerator.h"
#include "../../aux/Timer.h"

extern "C" {
#include "stinger.h"
}

namespace EnsembleClustering {

class GraphBenchmarkGTest: public testing::Test {

protected:

	int64_t n;

public:
	GraphBenchmarkGTest();
	virtual ~GraphBenchmarkGTest();
};

} /* namespace EnsembleClustering */
#endif /* GRAPHBENCHMARKGTEST_H_ */
