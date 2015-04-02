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

#include "../LinkPredictor.h"
#include "../../graph/Graph.h"

namespace NetworKit {

class LinkPredictionGTest : public testing::Test {

protected:
  Graph G;

  Graph trainingGraph;

  std::vector<std::pair<node, node>> missingLinks;

  std::vector<LinkPredictor::node_dyad_score_pair> predictions;

public:
  LinkPredictionGTest();

  void SetUp();

  //void TearDown();

};

} // namespace NetworKit

#endif /* LINKPREDICTIONGTEST_H_ */

#endif /* NOGTEST */