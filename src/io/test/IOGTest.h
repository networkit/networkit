/*
 * IOGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef IOGTEST_H_
#define IOGTEST_H_

#include <gtest/gtest.h>
#include <fstream>

#include "../../auxiliary/Log.h"
#include "../METISGraphReader.h"
#include "../METISGraphWriter.h"
#include "../ClusteringWriter.h"
#include "../ClusteringReader.h"
#include "../GraphIO.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../DotGraphWriter.h"

namespace NetworKit {

class IOGTest: public testing::Test {

};



} /* namespace NetworKit */
#endif /* IOGTEST_H_ */
