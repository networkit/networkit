/*
 * DynamicCommunityDetectionGTest.h
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef DYNAMICCOMMUNITYDETECTIONGTEST_H_
#define DYNAMICCOMMUNITYDETECTIONGTEST_H_

#include <gtest/gtest.h>

#include "../DynamicLabelPropagation.h"
#include "../../community/LabelPropagation.h"
#include "../../community/Louvain.h"
#include "../../generators/DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {

class DynamicCommunityDetectionGTest: public testing::Test {
public:
	DynamicCommunityDetectionGTest();
	virtual ~DynamicCommunityDetectionGTest();
};

} /* namespace NetworKit */
#endif /* DYNAMICCOMMUNITYDETECTIONGTEST_H_ */

#endif /*NOGTEST */
