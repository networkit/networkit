/*
 * CurveballGTest.h
 *
 *  Created on: Jul 13, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#ifndef RANDOMIZATION_TEST_CURVEBALL_H
#define RANDOMIZATION_TEST_CURVEBALL_H

#include <gtest/gtest.h>
#include "../../graph/Graph.h"

namespace NetworKit {

class CurveballGTest : public testing::Test  {
public:
    CurveballGTest() = default;
    virtual ~CurveballGTest() = default;

protected:
    void checkWithGraph(NetworKit::Graph&, bool checkBuilder = false);
};

}

#endif // RANDOMIZATION_TEST_CURVEBALL_H
