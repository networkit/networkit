/*
 * LocalDegreeGTest.h
 *
 *  Created on: 24.03.2015
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#ifndef LOCALDEGREETEST_H_
#define LOCALDEGREETEST_H_

#include <gtest/gtest.h>
#include "../../Globals.h"
#include "../../graph/Graph.h"

namespace NetworKit {

class LocalDegreeGTest: public testing::Test {
protected:
    static double getScore(const Graph& g, node x, node y, count rankX, count rankY);
};


} /* namespace NetworKit */
#endif /* LOCALDEGREETEST_H_ */

#endif /*NOGTEST */
