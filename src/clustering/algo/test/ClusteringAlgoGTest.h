/*
 * ClusteringAlgoGTest.h
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#ifndef CLUSTERINGALGOGTEST_H_
#define CLUSTERINGALGOGTEST_H_

#include <gtest/gtest.h>

#include "../../../aux/Log.h"
#include "../LabelPropagation.h"
#include "../../base/Modularity.h"
#include "../../../graph/GraphGenerator.h"
#include "../../../clustering/base/ClusteringGenerator.h"

namespace EnsembleClustering {

class ClusteringAlgoGTest: public testing::Test {
};

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGALGOGTEST_H_ */
