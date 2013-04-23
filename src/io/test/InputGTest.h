/*
 * InputGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef INPUTGTEST_H_
#define INPUTGTEST_H_

#include <gtest/gtest.h>
#include <fstream>

#include "../../auxiliary/Log.h"
#include "../METISGraphReader.h"
#include "../ClusteringWriter.h"
#include "../ClusteringReader.h"

#include "../GraphIO.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"

namespace NetworKit {

class InputGTest: public testing::Test {

};



} /* namespace NetworKit */
#endif /* INPUTGTEST_H_ */
