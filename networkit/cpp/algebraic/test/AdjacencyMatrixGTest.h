/*
 * AdjacencyMatrixGTest.h
 *
 *  Created on: 02.04.2014
 *      Author: Michael
 */

#ifndef NOGTEST

#ifndef ADJACENCYMATRIXGTEST_H_
#define ADJACENCYMATRIXGTEST_H_

#include "gtest/gtest.h"
#include "../AdjacencyMatrix.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"


namespace NetworKit {

class AdjacencyMatrixGTest : public testing::Test {
public:
	AdjacencyMatrixGTest();
	virtual ~AdjacencyMatrixGTest();
};


} /* namespace NetworKit */

#endif /* ADJACENCYMATRIXGTEST_H_ */

#endif
