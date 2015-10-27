/*
 * LCAGTest.cpp
 *
 *  Created on: 25.07.2014
 *      Author: dhoske
 */

#ifndef NOGTEST
/** @internal */

#include <unordered_set>
#include <random>

#include "LCAGTest.h"
#include "../Config.h"
#include "../RandomTreeGenerator.h"

namespace NetworKit {
namespace SDD {

/* Trivial LCA computation for comparison */
static node trivial_lca(const RootedTree& T, node u, node v) {
  std::unordered_set<node> nodes;
  node lca = node(-1);

  T.forNodesUpUntil([&] (node w) {
    nodes.emplace(w);
  }, u, T.getRoot());

  T.forNodesUpUntil([&] (node w) {
    if (lca == node(-1) && nodes.find(w) != nodes.end()) {
      lca = w;
    }
  }, v, T.getRoot());

  return lca;
}

/* Test on small manual tree */
TEST_F(LCAGTest, testManual) {
  /*
   *     0
   *    / \
   *   1   2
   *  / \   \
   * 5   3   4
   */
  RootedTree T({{0, 1, 1}, {1, 5, 1}, {1, 3, 1}, {0, 2, 1}, {2, 4, 1}}, 0);
  LCA lca(T);
  EXPECT_EQ(node(1), lca.query(5, 3));
  EXPECT_EQ(node(0), lca.query(4, 5));
  EXPECT_EQ(node(0), lca.query(1, 4));
  EXPECT_EQ(node(2), lca.query(4, 2));
  EXPECT_EQ(node(5), lca.query(5, 5));
  EXPECT_EQ(node(1), trivial_lca(T, 5, 3));
  EXPECT_EQ(node(0), trivial_lca(T, 4, 5));
  EXPECT_EQ(node(0), trivial_lca(T, 1, 4));
  EXPECT_EQ(node(2), trivial_lca(T, 4, 2));
  EXPECT_EQ(node(5), trivial_lca(T, 5, 5));
}

/* Test on large random graph. */
TEST_F(LCAGTest, testLarge) {
  static const SeedType SEED = 124234234;
  static const count N = (1 << 16);
  static const count TEST = 1000;

  RandomTreeGenerator gen(SEED);
  RootedTree T = gen.makeTree(N);
  LCA lca(T);

  // Test on random pairs of nodes
  std::default_random_engine engine(SEED);
  std::uniform_int_distribution<node> rand_node(0, N - 1);
  for (count i = 0; i < TEST; ++i) {
    node u = rand_node(engine), v = rand_node(engine);
    EXPECT_EQ(trivial_lca(T, u, v), lca.query(u, v));
  }
}

}
}

/** @endinternal */
#endif
