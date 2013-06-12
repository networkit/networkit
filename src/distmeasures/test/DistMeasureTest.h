/*
 * DistMeasureTest.h
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#ifndef DISTMEASURETEST_H_
#define DISTMEASURETEST_H_

#include <gtest/gtest.h>

#include "../AlgebraicDistances.h"
#include "../../graph/Graph.h"
#include "../../viz/PostscriptWriter.h"
#include "../../io/METISGraphReader.h"
#include "../../io/DibapGraphReader.h"
#include "../../clustering/Clustering.h"
#include "../../clustering/Modularity.h"
#include <cstdio>

namespace NetworKit {

class DistMeasureTest: public testing::Test {
public:
	DistMeasureTest();
	virtual ~DistMeasureTest();
};

} /* namespace NetworKit */
#endif /* DISTMEASURETEST_H_ */
