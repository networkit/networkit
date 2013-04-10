/*
 * GeneratorsTest.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef GENERATORSTEST_H_
#define GENERATORSTEST_H_

#include <gtest/gtest.h>

#include "../DynamicGraphGenerator.h"
#include "../DynamicBarabasiAlbertGenerator.h"
#include "../PubWebGenerator.h"
#include "../../viz/PostscriptWriter.h"
#include "../../clustering/ClusteringGenerator.h"

namespace NetworKit {

class GeneratorsTest: public testing::Test {
public:
	GeneratorsTest();
	virtual ~GeneratorsTest();
};

} /* namespace NetworKit */
#endif /* GENERATORSTEST_H_ */
