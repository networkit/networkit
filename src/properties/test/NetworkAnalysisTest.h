/*
 * NetworkAnalysisTest.h
 *
 *  Created on: 18.11.2013
 *      Author: christianocker
 */

#ifndef NOGTEST

#ifndef NETWORKANALYSISTEST_H_
#define NETWORKANALYSISTEST_H_

#include <gtest/gtest.h>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

#include "../GlobalClusteringCoefficient.h"
#include "../CoreDecomposition.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class NetworkAnalysisTest: public testing::Test {
public:
	NetworkAnalysisTest() {}
	virtual ~NetworkAnalysisTest() {}
};

} /* namespace NetworKit */
#endif /* NETWORKANALYSISTEST_H_ */

#endif /* NOGTEST */
