/*
 * DCDGTest.h
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef DCDGTEST_H_
#define DCDGTEST_H_

#include <gtest/gtest.h>

#include "../DynamicLabelPropagation.h"
#include "../../community/LabelPropagation.h"
#include "../../community/Louvain.h"
#include "../../generators/DynamicBarabasiAlbertGenerator.h"
#include "../../generators/DynamicPubWebGenerator.h"
#include "../../generators/DynamicDGSParser.h"
#include "../DynCDSetup.h"
#include "../PseudoDynamic.h"
#include "../../io/METISGraphReader.h"


namespace NetworKit {

class DCDGTest: public testing::Test {
public:
	DCDGTest();
	virtual ~DCDGTest();
};

} /* namespace NetworKit */
#endif /* DCDGTEST_H_ */

#endif /*NOGTEST */
