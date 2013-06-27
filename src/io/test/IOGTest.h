/*
 * IOGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef IOGTEST_H_
#define IOGTEST_H_

#include <gtest/gtest.h>
#include <fstream>
#include <unordered_set>
#include <vector>

#include "../../auxiliary/Log.h"
#include "../METISGraphReader.h"
#include "../METISGraphWriter.h"
#include "../ClusteringWriter.h"
#include "../ClusteringReader.h"
#include "../GraphIO.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../DotGraphWriter.h"
#include "../DGSReader.h"
#include "../EdgeListReader.h"
#include "../EdgeListClusteringReader.h"
#include "../SNAPEdgeListClusteringReader.h"
#include "../../clustering/Clustering.h"
#include "../../clustering/Modularity.h"
#include "../../community/LabelPropagation.h"



namespace NetworKit {

class IOGTest: public testing::Test {

};



} /* namespace NetworKit */
#endif /* IOGTEST_H_ */


#endif /* NOGTEST */
