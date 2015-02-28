/*
 * LinkPredictionGTest.cpp
 *
 *  Created on: 07.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NOGTEST

#include "LinkPredictionGTest.h"
#include "../../io/METISGraphReader.h"

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
  RandomEdgeRemover remover(G);
  std::pair<Graph, Graph> result = remover.remove(0.3);
  EXPECT_EQ(7, result.first.numberOfEdges());
  EXPECT_EQ(3, result.second.numberOfEdges());
}

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