/*
 * CycleDistributionTest.cpp
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef NOGTEST
/** @internal */

#include <memory>

#include "../../graph/GraphGenerator.h"
#include "../Config.h"
#include "CycleDistributionGTest.h"
using namespace std;

namespace NetworKit {
namespace SDD {

// Test settings
static const SeedType SEED = 1234;

/* Correctness test on an unweighted graph: complete graph with path as spanning tree. */
TEST_F(CycleDistributionGTest, testUniform) {
  GraphGenerator gen;
  static constexpr count SIZE = 100;

  Graph G = gen.makeCompleteGraph(SIZE);
  Graph GT = gen.makeCircularGraph(SIZE);
  GT.removeEdge(SIZE - 1, 0);
  RootedTree T(GT, 0);

  /* Off-tree-edges */
  vector<Edge> off_tree_edges;
  vector<edgeweight> stretches;
  for (index i = 0; i < SIZE; ++i) {
    for (index j = i + 2; j < SIZE; ++j) {
      off_tree_edges.emplace_back(i, j);
      stretches.emplace_back(1.0);
    }
  }

  for (int i = 1; i <= 2; ++i) {
    unique_ptr<CycleDistribution> dist;
    if (i == 1) {
      dist.reset(new UniformCycleDistribution(G, T, off_tree_edges, stretches, SEED));

      T.forEdges([&] (node u, node v) {
        EXPECT_EQ(0.0, dist->getProb({u, T.getParent(u)}));
      });
    } else {
      dist.reset(new StretchCycleDistribution(G, T, off_tree_edges, stretches, SEED));
    }

    for (int i = 0; i < 100; ++i) {
      Edge e = off_tree_edges[dist->randomCycle()];
      EXPECT_TRUE(G.hasEdge(e.u, e.v));
      EXPECT_FALSE(T.hasEdge(e.u, e.v));
    }
  }
}

/* Test of stretch probabilities on a small manual graph. */
TEST_F(CycleDistributionGTest, testStretchManual) {
  /* Not very pretty picture:
   *       0-------(2)
   *      / \        \
   *    (2) (3)      *3
   *    / ****\**(3)*
   *   1*´     2
   *   *      / \
   *  (1)   (1) (2)
   *   *    /     \
   *   `***4**(2)**5
   */
  Graph G = {
    /* ST */       {0, 1, 2}, {0, 2, 3}, {0, 3, 2}, {2, 4, 1}, {2, 5, 2},
    /* Off-tree */ {1, 3, 3}, {1, 4, 1}, {4, 5, 2}
  };
  RootedTree T({
    /* ST */       {0, 1, 2}, {0, 2, 3}, {0, 3, 2}, {2, 4, 1}, {2, 5, 2}
  }, 0);
  vector<pair<Edge, double>> expected_probs = {
    /* ST */       {{0, 1}, 0.0}, {{0, 2}, 0.0}, {{0, 3}, 0.0}, {{2, 4}, 0.0}, {{2, 5}, 0.0},
    /* Off-tree */ {{1, 3}, 3.0}, {{1, 4}, 11.0/6.0}, {{4, 5}, 3.0}
  };
  vector<Edge> off_tree_edges = {
    /* Off-tree */ {1, 3}, {1, 4}, {4, 5}
  };
  vector<edgeweight> stretches = {
    3.0, 11./6., 3.
  };

  StretchCycleDistribution dist(G, T, off_tree_edges, stretches, SEED);
  for (const auto& elem: expected_probs) {
    EXPECT_NEAR(elem.second, dist.getProb(elem.first), EPSILON);
  }
}

}
}

/** @endinternal */
#endif
