/*
 * CommonNeighboursIndex.h
 *
 *  Created on: 07.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NOGTEST

#ifndef LINKPREDICTIONGTEST_H_
#define LINKPREDICTIONGTEST_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../RandomEdgeRemover.h"

namespace NetworKit {

class LinkPredictionGTest : public testing::Test {

protected:
  Graph G;

public:
  LinkPredictionGTest();

  void SetUp();

  //void TearDown();

};

} // namespace NetworKit

#endif /* LINKPREDICTIONGTEST_H_ */

#endif /* NOGTEST */