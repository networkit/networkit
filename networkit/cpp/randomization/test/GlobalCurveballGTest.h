/*
 * GlobalCurveballGTest.h
 *
 *  Created on: 24.05.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef RANDOMIZATION_TEST_GLOBAL_CURVEBALL_GTEST_H
#define RANDOMIZATION_TEST_GLOBAL_CURVEBALL_GTEST_H

#include <gtest/gtest.h>
#include "../../graph/Graph.h"

namespace NetworKit {

class GlobalCurveballGTest : public ::testing::Test  {
public:
    GlobalCurveballGTest() = default;
    virtual ~GlobalCurveballGTest() = default;

protected:
    void checkWithGraph(NetworKit::Graph&);
};

}

#endif // RANDOMIZATION_TEST_GLOBAL_CURVEBALL_GTEST_H
