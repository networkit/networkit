/*
 * ClusteringAlgoGTest.h
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGALGOGTEST_H_
#define CLUSTERINGALGOGTEST_H_

#include <gtest/gtest.h>

#include "../../auxiliary/Log.h"
#include "../LabelPropagation.h"
#include "../Louvain.h"
#include "../../clustering/Modularity.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"

namespace NetworKit {

class ClusteringAlgoGTest: public testing::Test {
};

} /* namespace NetworKit */
#endif /* CLUSTERINGALGOGTEST_H_ */
