/*
 * ClusteringTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGGTEST_H_
#define CLUSTERINGGTEST_H_

#include <gtest/gtest.h>


#include "../../../aux/Log.h"
#include "../Clustering.h"
#include "../Modularity.h"
#include "../ModularitySequential.h"
#include "../Coverage.h"
#include "../ClusteringGenerator.h"
#include "../JaccardMeasure.h"
#include "../RandMeasure.h"
#include "../../../graph/GraphGenerator.h"
#include "../../algo/LabelPropagation.h"
#include "../../../io/METISGraphReader.h"
#include "../../../io/ClusteringReader.h"

namespace EnsembleClustering {

class ClusteringGTest: public testing::Test {


};



} /* namespace EnsembleClustering */
#endif /* CLUSTERINGTEST_H_ */
