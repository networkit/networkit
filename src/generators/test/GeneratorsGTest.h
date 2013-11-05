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

#include "../DynamicGraphSource.h"
#include "../DynamicBarabasiAlbertGenerator.h"
#include "../PubWebGenerator.h"
#include "../DynamicPubWebGenerator.h"
#include "../../viz/PostscriptWriter.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../community/LabelPropagation.h"
#include "../BTERGenerator.h"
#include "../../io/METISGraphWriter.h"
#include "../../io/DotGraphWriter.h"
#include "../BarabasiAlbertGenerator.h"
#include "../../io/GraphIO.h"
#include "../../io/METISGraphReader.h"
#include "../../properties/GraphProperties.h"

namespace NetworKit {

class GeneratorsGTest: public testing::Test {
public:
	GeneratorsGTest();
	virtual ~GeneratorsGTest();
};

} /* namespace NetworKit */
#endif /* GENERATORSGTEST_H_ */

#endif /*NOGTEST */
