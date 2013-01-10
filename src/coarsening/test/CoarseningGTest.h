/*
 * CoarseningGTest.h
 *
 *  Created on: 20.12.2012
 *      Author: cls
 */

#ifndef COARSENINGGTEST_H_
#define COARSENINGGTEST_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../coarsening/ClusterContracter.h"
#include "../../coarsening/GraphContraction.h"

namespace EnsembleClustering {

/**
 * googletest test fixture for the coarsening module.
 */
class CoarseningGTest: public testing::Test {

	// TODO: are constructor/destructor needed?

};


} /* namespace EnsembleClustering */
#endif /* COARSENINGGTEST_H_ */
