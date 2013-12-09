/*
 * ConnectedComponentsGTest.h
 *
 *  Created on: Sep 16, 2013
 *      Author: birderon
 */

#ifndef NOGTEST

#ifndef CONNECTEDCOMPONENTSGTEST_H_
#define CONNECTEDCOMPONENTSGTEST_H_

#include <gtest/gtest.h>
#include "../ConnectedComponents.h"
//#include "../../io/METISGraphReader.h"

namespace NetworKit {

class ConnectedComponentsGTest: public testing::Test {
public:
	ConnectedComponentsGTest();
	virtual ~ConnectedComponentsGTest();
};

} /* namespace NetworKit */
#endif /* CONNECTEDCOMPONENTSGTEST_H_ */

#endif /*NOGTEST */







