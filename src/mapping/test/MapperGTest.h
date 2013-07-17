/*
 * MapperGTest.h
 *
 *  Created on: Jul 15, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#ifndef MAPPERGTEST_H_
#define MAPPERGTEST_H_

#include <gtest/gtest.h>

#include "../StaticMapper.h"
#include "../RecBisMapper.h"
#include "../RcmMapper.h"
#include "../RcmMapperWW.h"
#include "../GreedyMapper.h"
#include "../../partitioning/BalancedLabelPropagation.h"
#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../io/DibapGraphReader.h"
#include "../../io/METISGraphReader.h"
#include "../../io/METISGraphWriter.h"
#include "../../clustering/Clustering.h"
#include "../../clustering/ClusteringGenerator.h"


namespace NetworKit {

class MapperGTest: public testing::Test {
public:
	MapperGTest();
	virtual ~MapperGTest();
};

}

#endif /* MAPPERGTEST_H_ */

#endif
