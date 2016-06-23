/*
 * Graph2Benchmark.h
 *
 *  Created on: 05.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef GRAPH2BENCHMARK_H_
#define GRAPH2BENCHMARK_H_

#include <gtest/gtest.h>

#include "../Graph.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Log.h"

namespace NetworKit {

class Graph2Benchmark: public testing::Test {
public:
	Graph2Benchmark();
	virtual ~Graph2Benchmark();
};

} /* namespace NetworKit */
#endif /* GRAPH2BENCHMARK_H_ */

#endif /*NOGTEST */
