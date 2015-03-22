/*
 * LinkPredictionGTest.cpp
 *
 *  Created on: 07.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NOGTEST

#include "LinkPredictionGTest.h"
#include "../../io/METISGraphReader.h"
#include "../KatzIndex.h"
#include "../CommonNeighborsIndex.h"
#include "../EdgeSelector.h"
#include "../ROCMetric.h"
#include "../PrecisionRecallMetric.h"
#include "../RandomEdgePartitioner.h"
#include "../KFoldCrossValidator.h"
#include "../UnconnectedNodesFinder.h"

namespace NetworKit {

LinkPredictionGTest::LinkPredictionGTest() : G(7) {
}

void LinkPredictionGTest::SetUp() {
  G.addEdge(0, 1);
  G.addEdge(0, 3);
  G.addEdge(1, 2);
  G.addEdge(1, 4);
  G.addEdge(2, 3);
  G.addEdge(2, 4);
  G.addEdge(2, 5);
  G.addEdge(3, 4);
  G.addEdge(3, 5);
  G.addEdge(4, 5);
}

TEST_F(LinkPredictionGTest, testRandomEdgeRemoval) {
  RandomEdgePartitioner partitioner(G);
  std::pair<Graph, Graph> graphPartitions = partitioner.partitionByPercentage(0.3);

  EXPECT_EQ(7, graphPartitions.first.numberOfEdges());
  EXPECT_EQ(3, graphPartitions.second.numberOfEdges());
}

TEST_F(LinkPredictionGTest, testEdgeSelectorGetByCount) {
  LinkPredictor* predictor = new KatzIndex(G, 2, 1);
  EdgeSelector selector(G, predictor);
  selector.calculateScores();
  std::vector<std::pair<std::pair<node, node>, double>> edges = selector.getByCount(6);
  // Expects that the edges are ordered descendingly by scores and that on equality the edge-pairs
  // are ordered ascendingly [(0,1) < (0,2) and (1, 2) < (1, 0)].
  EXPECT_EQ(1, edges[0].first.first); EXPECT_EQ(3, edges[0].first.second); EXPECT_EQ(3, edges[0].second);
  EXPECT_EQ(0, edges[1].first.first); EXPECT_EQ(2, edges[1].first.second); EXPECT_EQ(2, edges[1].second);
  EXPECT_EQ(0, edges[2].first.first); EXPECT_EQ(4, edges[2].first.second); EXPECT_EQ(2, edges[2].second);
  EXPECT_EQ(1, edges[3].first.first); EXPECT_EQ(5, edges[3].first.second); EXPECT_EQ(2, edges[3].second);
  EXPECT_EQ(0, edges[4].first.first); EXPECT_EQ(5, edges[4].first.second); EXPECT_EQ(1, edges[4].second);
  EXPECT_EQ(0, edges[5].first.first); EXPECT_EQ(6, edges[5].first.second); EXPECT_EQ(0, edges[5].second);
  delete predictor;
}

TEST_F(LinkPredictionGTest, testEdgeSelectorGetByCountArgumentException) {
  LinkPredictor* predictor = new KatzIndex(G, 2, 1);
  EdgeSelector selector(G, predictor, 5);
  selector.calculateScores();
  EXPECT_THROW(selector.getByCount(6), std::invalid_argument);
  delete predictor;
}

TEST_F(LinkPredictionGTest, testEdgeSelectorGetAllMissingCalcCall) {
  LinkPredictor* predictor = new KatzIndex(G, 2, 1);
  EdgeSelector selector(G, predictor, 5);
  EXPECT_THROW(selector.getAll(), std::logic_error);
  delete predictor;
}

TEST_F(LinkPredictionGTest, testEdgeSelectorGetByCountMissingCalcCall) {
  LinkPredictor* predictor = new KatzIndex(G, 2, 1);
  EdgeSelector selector(G, predictor);
  EXPECT_THROW(selector.getByCount(1), std::logic_error);
  delete predictor;
}

/*TEST_F(LinkPredictionGTest, testKFoldCrossValidator) {
  KatzIndex katzIndex;
  METISGraphReader graphReader;
  Graph jazz = graphReader.read("input/jazz.graph");
  ROC roc;

  KFoldCrossValidator validator(jazz, &katzIndex, &roc);
  double average = validator.crossValidate(10);
  EXPECT_NEAR(average, 0.78, 0.03);
}*/

/*TEST_F(LinkPredictionGTest, testMissingLinkFinder) {
  METISGraphReader graphReader;
  Graph jazz = graphReader.read("input/PGPgiantcompo.graph");
  UnconnectedNodesFinder unf(jazz);
  std::vector<std::pair<node, node>> missingEdges = unf.findAll(2);
  INFO(missingEdges.size());
}*/

TEST_F(LinkPredictionGTest, testPrecisionRecallMetric) {
  RandomEdgePartitioner partitioner(G);
  std::pair<Graph, Graph> graphPartitions = partitioner.partitionByPercentage(0.3);
  CommonNeighborsIndex cni(graphPartitions.first);
  std::vector<std::pair<std::pair<node, node>, double>> results = cni.runAll();

  INFO("#Predictions: ", results.size());

  PrecisionRecallMetric pr(graphPartitions.second, results);
  pr.generatePoints();
  EXPECT_NEAR(0.7, pr.areaUnderCurve(), 0.01);
}

/*TEST_F(LinkPredictionGTest, testReceiverOperatingCharacteristic) {
  RandomEdgePartitioner partitioner(G);
  std::pair<Graph, Graph> graphPartitions = partitioner.partition(0.3);

  KatzIndex katz(graphPartitions.first, 2, 1);
  std::vector<LinkPredictor::node_dyad_score_pair> scores = katz.runAll(4);
  for (index i = 0; i < scores.size(); ++i) {
    INFO("entries[", i, "] = ((", scores[i].first.first, ", ", scores[i].first.second, "), ", scores[i].second, ")");
  }
}*/


/*
TEST_F(LinkPredictionGTest, testCommonNeighborsRun) {
  CommonNeighborsIndex cni(G);
  EXPECT_EQ(3.0, cni.run(1, 3));
  EXPECT_EQ(2.0, cni.run(3, 5));
  EXPECT_EQ(0.0, cni.run(4, 6));
}

TEST_F(LinkPredictionGTest, testCommonNeighborsRunAll) {
  CommonNeighborsIndex cni(G);
  ScoreCollection scores = cni.runAll();
  EXPECT_EQ(3.0, scores.getScore(1, 3));
  EXPECT_EQ(2.0, scores.getScore(3, 5));
  EXPECT_EQ(0.0, scores.getScore(4, 6));
}*/
/*
TEST_F(LinkPredictionGTest, testPreferentialAttachmentRun) {
  PreferentialAttachmentIndex pai(G);
  EXPECT_EQ(12.0, pai.run(1, 3));
  EXPECT_EQ(16.0, pai.run(2, 3));
  EXPECT_EQ(0.0, pai.run(4, 6));
}

TEST_F(LinkPredictionGTest, testPreferentialAttachmentRunAll) {
  PreferentialAttachmentIndex pai(G);
  ScoreCollection scores = pai.runAll();
  EXPECT_EQ(12.0, scores.getScore(1, 3));
  EXPECT_EQ(16.0, scores.getScore(2, 3));
  EXPECT_EQ(0.0, scores.getScore(4, 6));
}

TEST_F(LinkPredictionGTest, testJaccardCoefficientRun) {
  JaccardCoefficientIndex jci(G);
  EXPECT_DOUBLE_EQ(0.75, jci.run(1, 3));
  EXPECT_DOUBLE_EQ((double) 1 / 3, jci.run(2, 3));
  EXPECT_EQ(0.0, jci.run(4, 6));
}

TEST_F(LinkPredictionGTest, testJaccardCoefficientRunAll) {
  JaccardCoefficientIndex jci(G);
  ScoreCollection scores = jci.runAll();
  EXPECT_DOUBLE_EQ(0.75, scores.getScore(1, 3));
  EXPECT_DOUBLE_EQ((double) 1 / 3, scores.getScore(2, 3));
  EXPECT_EQ(0.0, scores.getScore(4, 6));
}
*/
/*TEST_F(LinkPredictionGTest, testNonexistentResultsException) {
  PreferentialAttachmentIndex pai(G);
  EXPECT_THROW(pai.getResult(1, 2), std::runtime_error);
}*/

// TODO: Write a test to make sure the score of a node with itself
// is [...] (maybe 0, MISSING_SCORE or a new INVALID_SCORE)
// Maybe it should throw an exception?

} // namespace NetworKit

#endif /* NOGTEST */