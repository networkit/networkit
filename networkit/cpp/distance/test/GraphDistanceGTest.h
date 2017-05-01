/*
 * GraphDistanceGTest.h
 *
 *  Created on: 23.07.2013
 *      Author: Henning Meyerhenke (meyerhenke@kit.edu)
 */

#ifndef NOGTEST

#ifndef GRAPHDISTANCEGTEST_H_
#define GRAPHDISTANCEGTEST_H_

#include <gtest/gtest.h>

#include "../GraphDistance.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class GraphDistanceGTest: public testing::Test {
public:
	GraphDistanceGTest();
	virtual ~GraphDistanceGTest();
};

} /* namespace NetworKit */
#endif /* GRAPHDISTANCEGTEST_H_ */

#endif /*NOGTEST */
