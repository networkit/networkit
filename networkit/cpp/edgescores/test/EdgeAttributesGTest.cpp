/*
 *
 */

#ifndef NOGTEST

#include "EdgeAttributesGTest.h"

#include "../../io/METISGraphReader.h"
#include "../PrefixJaccardCoefficient.h"
#include "../../sparsification/SimmelianJaccardAttributizer.h"
#include "../ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

TEST_F(EdgeAttributesGTest, testPrefixJaccardCoefficient) {
	Graph karate = METISGraphReader().read("input/karate.graph");
	karate.indexEdges();
	std::vector<count> triangles = ChibaNishizekiTriangleCounter(karate).getAttribute();

	PrefixJaccardCoefficient<count> prefixJaccard(karate, triangles);
	prefixJaccard.run();

	std::vector<double> simmelianJaccard = SimmelianJaccardAttributizer(karate, triangles).getAttribute();

	EXPECT_EQ(simmelianJaccard, prefixJaccard.getAttribute());
}

}

#endif