/*
 * PropertiesGTest.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef PROPERTIESGTEST_H_
#define PROPERTIESGTEST_H_

#include <gtest/gtest.h>

#include "../ClusteringCoefficient.h"
#include "../../graph/GraphGenerator.h"
#include "../../properties/GraphProperties.h"
#include "../../io/METISGraphReader.h"


namespace NetworKit {

class PropertiesGTest: public testing::Test {
public:
	PropertiesGTest();
	virtual ~PropertiesGTest();
};

} /* namespace NetworKit */
#endif /* PROPERTIESGTEST_H_ */
