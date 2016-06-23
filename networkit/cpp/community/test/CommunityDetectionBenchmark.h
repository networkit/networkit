/*
 * CommunityDetectionBenchmark.h
 *
 *  Created on: 16.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef COMMUNITY_DETECTION_BENCHMARK_H_
#define COMMUNITY_DETECTION_BENCHMARK_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class CommunityDetectionBenchmark: public testing::Test {
public:
	virtual void SetUp();

protected:
	METISGraphReader metisReader;	

};

} /* namespace NetworKit */
#endif /* COMMUNITY_DETECTION_BENCHMARK_H_ */

#endif /*NOGTEST */
