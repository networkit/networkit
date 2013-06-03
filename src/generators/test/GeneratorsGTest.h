/*
 * GeneratorsGTest.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef GENERATORSGTEST_H_
#define GENERATORSGTEST_H_

#include <gtest/gtest.h>

#include "../DynamicGraphGenerator.h"
#include "../DynamicBarabasiAlbertGenerator.h"
#include "../PubWebGenerator.h"
#include "../../viz/PostscriptWriter.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../community/LabelPropagation.h"
#include "../BTERGenerator.h"
#include "../../io/METISGraphWriter.h"
#include "../StaticBarabasiAlbertGenerator.h"
#include "../../io/GraphIO.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class GeneratorsGTest: public testing::Test {
public:
	GeneratorsGTest();
	virtual ~GeneratorsGTest();
};

} /* namespace NetworKit */
#endif /* GENERATORSGTEST_H_ */
