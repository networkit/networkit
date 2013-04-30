/*
 * GraphBenchmark.h
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHBENCHMARK_H_
#define GRAPHBENCHMARK_H_

#include <gtest/gtest.h>

#include "../GraphGenerator.h"
#include "../../auxiliary/Timer.h"


namespace NetworKit {

class GraphBenchmark: public testing::Test {
protected:
	int64_t n;
public:
	GraphBenchmark();
	virtual ~GraphBenchmark();
};

} /* namespace NetworKit */
#endif /* GRAPHBENCHMARK_H_ */
