/*
 * BetweennessCentralityGTest.h
 *
 *  Created on: 09.12.2013
 *      Author: Lukas Barth, David Weiß
 */

#ifndef NOGTEST

#ifndef BetweennessCentralityGTest_H_
#define BetweennessCentralityGTest_H_

#include <gtest/gtest.h>
#include "../BetweennessCentrality.h"
//#include "../../io/METISGraphReader.h"

namespace GrauBart {

class BetweennessCentralityGTest: public testing::Test {
public:
	BetweennessCentralityGTest();
	virtual ~BetweennessCentralityGTest();
};

} /* namespace GrauBart */
#endif /* BetweennessCentralityGTest_H_ */

#endif /*NOGTEST */







