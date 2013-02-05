/*
 * Graph2Benchmark.h
 *
 *  Created on: 05.02.2013
 *      Author: cls
 */

#ifndef GRAPH2BENCHMARK_H_
#define GRAPH2BENCHMARK_H_

#include <gtest/gtest.h>

#include "../Graph2.h"
#include "../../aux/Timer.h"
#include "../../aux/Log.h"

namespace EnsembleClustering {

class Graph2Benchmark: public testing::Test {
public:
	Graph2Benchmark();
	virtual ~Graph2Benchmark();
};

} /* namespace EnsembleClustering */
#endif /* GRAPH2BENCHMARK_H_ */
