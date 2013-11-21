/*
 * PropertiesGTest_Ritter.h
 */

#ifndef NOGTEST

#ifndef PROPERTIESGTEST_RITTER_H_
#define PROPERTIESGTEST_RITTER_H_

#include <gtest/gtest.h>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

#include "../ClusteringCoefficient.h"
#include "../../graph/GraphGenerator.h"
#include "../../properties/GraphProperties.h"
#include "../../io/METISGraphReader.h"


namespace NetworKit {

class PropertiesGTest_Ritter: public testing::Test {
public:
	PropertiesGTest_Ritter();
	virtual ~PropertiesGTest_Ritter();
};

} /* namespace NetworKit */
#endif /* PROPERTIESGTEST_RITTER_H_ */

#endif /* NOGTEST */
