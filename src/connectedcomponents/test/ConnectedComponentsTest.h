/*
 * findCCGTest.h
 *
 *  Created on: Sep 16, 2013
 *      Author: birderon
 */

#ifndef NOGTEST

#ifndef CONNECTEDCOMPONENTSTEST_H_
#define CONNECTEDCOMPONENTSTEST_H_

#include <gtest/gtest.h>
#include "../ConnectedComponents.h"
//#include "../../io/METISGraphReader.h"

namespace NetworKit {

class ConnectedComponentsTest: public testing::Test {
public:
	ConnectedComponentsTest();
	virtual ~ConnectedComponentsTest();
};

} /* namespace NetworKit */
#endif /* CONNECTEDCOMPONENTSTEST_H_ */

#endif /*NOGTEST */







