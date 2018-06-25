/*
 * GlobalCurveballBenchmark.h
 *
 *  Created on: 24.05.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef RANDOMIZATION_TEST_GLOBAL_CURVEBALL_BENCHMARK_H
#define RANDOMIZATION_TEST_GLOBAL_CURVEBALL_BENCHMARK_H

#include <gtest/gtest.h>
#include "../../graph/Graph.h"

namespace NetworKit {

class GlobalCurveballBenchmark : public ::testing::Test  {
public:
    GlobalCurveballBenchmark() = default;
    virtual ~GlobalCurveballBenchmark() = default;

protected:
    void checkWithGraph(NetworKit::Graph&);
};

}

#endif // RANDOMIZATION_TEST_GLOBAL_CURVEBALL_BENCHMARK_H
