/*
 * CoarseningGTest.h
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef COARSENINGGTEST_H_
#define COARSENINGGTEST_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/base/ClusteringGenerator.h"
#include "../../coarsening/ClusterContracter.h"
#include "../../coarsening/GraphContraction.h"
#include "../../coarsening/ClusteringProjector.h"

namespace EnsembleClustering {

/**
 * googletest test fixture for the coarsening module.
 */
class CoarseningGTest: public testing::Test {

	// TODO: are constructor/destructor needed?

};


} /* namespace EnsembleClustering */
#endif /* COARSENINGGTEST_H_ */
