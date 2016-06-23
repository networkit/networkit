/*
 * SpanningEdgeCentralityGTest.h
 *
 *  Created on: Jan 17, 2016
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_CENTRALITY_TEST_SPANNINGEDGECENTRALITYGTEST_H_
#define NETWORKIT_CPP_CENTRALITY_TEST_SPANNINGEDGECENTRALITYGTEST_H_

#include "gtest/gtest.h"
#include "../SpanningEdgeCentrality.h"

#include <vector>
#include <string>

namespace NetworKit {

using namespace std;

class SpanningEdgeCentralityGTest : public testing::Test {
public:
	SpanningEdgeCentralityGTest() = default;
	virtual ~SpanningEdgeCentralityGTest() = default;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_CENTRALITY_TEST_SPANNINGEDGECENTRALITYGTEST_H_ */
