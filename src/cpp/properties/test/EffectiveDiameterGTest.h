/*
 * EffectiveDiameterGTest.h
 *
 *  Created on: Jun 16, 2014
 *      Author: Marc Nemes
 */

#ifndef NOGTEST

#ifndef EFFECTIVEDIAMETERGTEST_H_
#define EFFECTIVEDIAMETERGTEST_H_

#include <gtest/gtest.h>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <fstream> // for ofstream
#include <string>
#include <list>

#include "../../graph/GraphGenerator.h"
#include "../../properties/GraphProperties.h"
#include "../../io/METISGraphReader.h"


namespace NetworKit {

class EffectiveDiameterGTest: public testing::Test {
public:
	EffectiveDiameterGTest();
	virtual ~EffectiveDiameterGTest();
};

} /* namespace NetworKit */
#endif /* EFFECTIVEDIAMETERGTEST_H_ */

#endif /* NOGTEST */
