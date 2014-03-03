/*
 * MatcherGTest.h
 *
 *  Created on: Jun 14, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#ifndef MATCHERGTEST_H_
#define MATCHERGTEST_H_

#include <gtest/gtest.h>

#include "../Matcher.h"
#include "../Matching.h"
#include "../PathGrowingMatcher.h"
#include "../ParallelMatcher.h"
#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../io/DibapGraphReader.h"


namespace NetworKit {

class MatcherGTest: public testing::Test {
public:
	MatcherGTest();
	virtual ~MatcherGTest();
};

}

#endif /* MATCHERGTEST_H_ */

#endif
