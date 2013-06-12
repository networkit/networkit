/*
 * EnsembleGTest.h
 *
 *  Created on: 31.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#ifndef ENSEMBLEGTEST_H_
#define ENSEMBLEGTEST_H_

#include <gtest/gtest.h>


#include "../EnsembleMultilevel.h"
#include "../EnsemblePreprocessing.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/Modularity.h"
#include "../../clustering/JaccardMeasure.h"
#include "../../clustering/NodeStructuralRandMeasure.h"
#include "../../community/LabelPropagation.h"
#include "../../community/Louvain.h"
#include "../../overlap/HashingOverlapper.h"


namespace NetworKit {

class EnsembleGTest: public testing::Test {

};




} /* namespace NetworKit */
#endif /* ENSEMBLEGTEST_H_ */

#endif /*NOGTEST */
