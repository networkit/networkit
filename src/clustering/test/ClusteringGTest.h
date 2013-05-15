/*
 * ClusteringTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGGTEST_H_
#define CLUSTERINGGTEST_H_

#include <gtest/gtest.h>


#include "../../auxiliary/Log.h"
#include "../Clustering.h"
#include "../Modularity.h"
#include "../ModularitySequential.h"
#include "../Coverage.h"
#include "../ClusteringGenerator.h"
#include "../JaccardMeasure.h"
#include "../NodeStructuralRandMeasure.h"
#include "../GraphStructuralRandMeasure.h"
#include "../../graph/GraphGenerator.h"
#include "../../io/METISGraphReader.h"
#include "../../io/ClusteringReader.h"
#include "../ClusteringCoefficient.h"

namespace NetworKit {

class ClusteringGTest: public testing::Test {


};



} /* namespace NetworKit */
#endif /* CLUSTERINGTEST_H_ */
