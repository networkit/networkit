/*
 * DistMeasureTest.h
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */
#ifndef NOGTEST

#ifndef DISTMEASURESGTEST_H_
#define DISTMEASURESGTEST_H_

#include <gtest/gtest.h>
#include <cstdio>


#include "../../graph/Graph.h"
#include "../../viz/PostscriptWriter.h"
#include "../../io/METISGraphReader.h"
#include "../../io/DibapGraphReader.h"
#include "../../structures/Partition.h"
#include "../../community/Modularity.h"
#include "../AlgebraicDistanceIndex.h"

namespace NetworKit {

class DistMeasuresGTest: public testing::Test {
public:
	DistMeasuresGTest();
	virtual ~DistMeasuresGTest();
};

} /* namespace NetworKit */
#endif /* DISTMEASURETEST_H_ */

#endif /*NOGTEST */
