/*
 * ClusteringAlgoGTest.h
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGALGOGTEST_H_
#define CLUSTERINGALGOGTEST_H_

#include <gtest/gtest.h>

#include "../../../aux/Log.h"
#include "../LabelPropagation.h"
#include "../Louvain.h"
#include "../../base/Modularity.h"
#include "../../../graph/GraphGenerator.h"
#include "../../../clustering/base/ClusteringGenerator.h"

namespace EnsembleClustering {

class ClusteringAlgoGTest: public testing::Test {
};

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGALGOGTEST_H_ */
