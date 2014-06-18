/*
 * EffectiveDiameterGTest.cpp
 *
 *  Created on: Jun 16, 2014
 *      Author: Marc Nemes
 */

#ifndef NOGTEST

#include "EffectiveDiameterGTest.h"
#include "../EffectiveDiameter.h"


namespace NetworKit {

EffectiveDiameterGTest::EffectiveDiameterGTest() {

}

EffectiveDiameterGTest::~EffectiveDiameterGTest() {

}
TEST_F(EffectiveDiameterGTest, testEffectiveDiameter) {
	using namespace std;
	vector<pair<string, int>> testInstances= {
		pair<string, int>("celegans_metabolic", 6),
		pair<string, int>("jazz", 5),
		pair<string, int>("karate", 5),
		pair<string, int>("lesmis", 5)
	};

	for (auto testInstance : testInstances) {
		METISGraphReader reader;
		Graph G = reader.read("input/" + testInstance.first + ".graph");
		int diameter = EffectiveDiameter::effectiveDiameter(G);
		EXPECT_EQ(diameter, testInstance.second);
	}
}
} /* namespace NetworKit */

#endif /*NOGTEST*/
