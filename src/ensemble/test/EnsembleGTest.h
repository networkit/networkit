/*
 * EnsembleGTest.h
 *
 *  Created on: 31.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ENSEMBLEGTEST_H_
#define ENSEMBLEGTEST_H_

#include <gtest/gtest.h>


#include "../EnsembleClusterer.h"
#include "../EnsemblePreprocessing.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/base/Modularity.h"
#include "../../clustering/base/JaccardMeasure.h"
#include "../../clustering/base/RandMeasure.h"
#include "../../clustering/algo/LabelPropagation.h"
#include "../../clustering/algo/Louvain.h";
#include "../../overlap/HashingOverlapper.h"


namespace EnsembleClustering {

class EnsembleGTest: public testing::Test {

};




} /* namespace EnsembleClustering */
#endif /* ENSEMBLEGTEST_H_ */
