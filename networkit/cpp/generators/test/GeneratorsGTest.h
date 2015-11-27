/*
 * GeneratorsGTest.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef GENERATORSGTEST_H_
#define GENERATORSGTEST_H_

#include <gtest/gtest.h>

#include "../HyperbolicGenerator.h"
#include "../DynamicHyperbolicGenerator.h"

namespace NetworKit {

class GeneratorsGTest: public testing::Test {
public:
	GeneratorsGTest();

	vector<double> getAngles(DynamicHyperbolicGenerator dynGen) {
		return dynGen.angles;
	}

	vector<double> getRadii(DynamicHyperbolicGenerator dynGen) {
		return dynGen.radii;
	}

	count getQuadTreeHeight(DynamicHyperbolicGenerator dynGen) {
		return dynGen.quadTreeHeight();
	}
};

} /* namespace NetworKit */
#endif /* GENERATORSGTEST_H_ */

#endif /*NOGTEST */
