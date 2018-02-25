/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "DynTopClosenessGTest.h"
#include "../../generators/DorogovtsevMendesGenerator.h"
#include "../DynTopHarmonicCloseness.h"

namespace NetworKit {

TEST_F(DynTopClosenessGTest, testDynTopHarmonicCloseness) {

  Aux::Random::setSeed(42, false);
  Graph G = DorogovtsevMendesGenerator(500).generate();

  count k = 10;

  DynTopHarmonicCloseness centrality(G, k, true);

  centrality.run();

  count numInsertions = 100;

  std::vector<GraphEvent> deletions;
  std::vector<GraphEvent> insertions;

  for (count i = 0; i < numInsertions; i++) {

    node u = G.upperNodeIdBound();
    node v = G.upperNodeIdBound();

    do {
      u = G.randomNode();
      v = G.randomNode();

    } while (G.hasEdge(u, v));

    GraphEvent edgeAddition(GraphEvent::EDGE_ADDITION, u, v);
    insertions.insert(insertions.begin(), edgeAddition);

    GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
    deletions.push_back(edgeDeletion);

    G.addEdge(u, v);
  }

  for (auto e : insertions) {
    G.removeEdge(e.u, e.v);
  }

  for (GraphEvent edgeAddition : insertions) {

    node u = edgeAddition.u;
    node v = edgeAddition.v;

    G.addEdge(u, v);
    centrality.update(edgeAddition);

    DynTopHarmonicCloseness reference(G, k, false);
    reference.run();

    auto scores = centrality.ranking();
    auto referenceScores = reference.ranking();

    for (count j = 0; j < k; j++) {
      bool nodeMatch = referenceScores[j].first == scores[j].first;
      bool scoreMatch =
          std::abs(scores[j].second - referenceScores[j].second) <= 0.01;
      if (!scoreMatch) {
        INFO(scores[j].second, " - ", referenceScores[j].second);
      }
      EXPECT_TRUE(nodeMatch && scoreMatch);
    }
  }

  std::reverse(deletions.begin(), deletions.end());

  for (auto deletion : deletions) {
    G.removeEdge(deletion.u, deletion.v);
    centrality.update(deletion);

    DynTopHarmonicCloseness reference(G, k);
    reference.run();

    auto scores = centrality.ranking();
    auto referenceScores = reference.ranking();

    for (count j = 0; j < k; j++) {
      bool nodeMatch = referenceScores[j].first == scores[j].first;
      bool scoreMatch =
          std::abs(scores[j].second - referenceScores[j].second) <= 0.01;

      EXPECT_TRUE(nodeMatch && scoreMatch);
    }
  }
}
} // namespace NetworKit
