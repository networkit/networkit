/*
 * EnsembleGTest.h
 *
 *  Created on: 31.12.2012
 *      Author: cls
 */

#ifndef ENSEMBLEGTEST_H_
#define ENSEMBLEGTEST_H_

#include <gtest/gtest.h>


#include "../EnsembleClusterer.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/Modularity.h"
#include "../../clustering/LabelPropagation.h"

namespace EnsembleClustering {

class EnsembleGTest: public testing::Test {

};




} /* namespace EnsembleClustering */
#endif /* ENSEMBLEGTEST_H_ */
