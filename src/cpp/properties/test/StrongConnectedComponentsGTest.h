/*
 * StrongConnectedComponentsGTest.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef STRONGCONNECTEDCOMPONENTSGTEST_H_
#define STRONGCONNECTEDCOMPONENTSGTEST_H_

#include <gtest/gtest.h>

#include "../../structures/Partition.h"

namespace NetworKit {

class StrongConnectedComponentsGTest: public testing::Test {
public:
	StrongConnectedComponentsGTest();
	virtual ~StrongConnectedComponentsGTest();

protected:
	void comparePartitions(const Partition& p1, const Partition& p2) const;
};

} /* namespace NetworKit */
#endif /* STRONGCONNECTEDCOMPONENTSGTEST_H_ */

#endif /*NOGTEST */







