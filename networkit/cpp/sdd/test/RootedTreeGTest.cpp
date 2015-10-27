/*
 * RootedTreeGTest.cpp
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef NOGTEST
/** @internal */

#include "RootedTreeGTest.h"
#include "../Config.h"
#include "../RootedTree.h"
#include "../RandomTreeGenerator.h"

using namespace std;

namespace NetworKit {
namespace SDD {

/* Is forNodesPost a post-order traversal? */
TEST_F(RootedTreeGTest, testTraverse) {
  RootedTree t(Graph {
    {1, 0, 0.1},
    {1, 2, 0.1},
    {2, 3, 0.1},
    {2, 4, 0.1}
  }, 1);

  vector<node> expected = {0, 3, 4, 2, 1};
  vector<node> actual;
  t.forNodesPost([&] (node u) {
    actual.emplace_back(u);
  });

  EXPECT_EQ(expected, actual);
}

/* Is upwards-traversal correct? */
TEST_F(RootedTreeGTest, testTraverseUp) {
  RootedTree t(Graph {
    {1, 0, 1},
    {1, 2, 2},
    {2, 3, 5},
    {2, 4, 3}
  }, 1);

  edgeweight expected = 5+2;
  edgeweight actual = 0;
  t.forEdgesUp([&] (node u, node v, edgeweight weight) {
    actual += weight;
  }, 3);

  EXPECT_EQ(expected, actual);
}

/* Is Eulerian traversal correct? */
TEST_F(RootedTreeGTest, testTraverseEulerian) {
  RootedTree t(Graph {
    {0, 1, 1},
    {1, 2, 1},
    {2, 3, 1},
  }, 0);
  std::vector<pair<node, count>> expected = {
     {0, 0}, {1, 1}, {2, 2}, {3, 3}, {2, 2}, {1, 1}, {0, 0}
  };

  index i = 0;
  t.forNodesEulerian([&] (node u, count depth) {
    ASSERT_LT(i, expected.size());
    ASSERT_EQ(expected[i++], make_pair(u, depth));
  });
  EXPECT_EQ(expected.size(), i);
}

/* Is the list of children set up correctly? */
TEST_F(RootedTreeGTest, testGetChildren) {
  RootedTree t(Graph {
    {1, 0, 1},
    {0, 5, 1},
    {0, 6, 2},
    {1, 2, 0},
    {2, 3, 4},
    {2, 4, 0}
  }, 1);

  edgeweight expected = 1+2+4+0;
  edgeweight actual = 0;
  t.forNodesPost([&] (node u) {
    if (u % 2 == 0) {
      for (const auto& v: t.getChildren(u)) {
        actual += t.getWeight(u, v);
      }
    }
  });

  EXPECT_EQ(expected, actual);
}

/* Are random trees generated correctly? */
TEST_F(RootedTreeGTest, testRandomTree) {
  static const count n = 1000;
  static const edgeweight lower = 1.0, upper = 2.0;
  RandomTreeGenerator gen(1234);

  auto T1 = gen.makeTree(n, lower, upper);
  EXPECT_EQ(n, T1.numberOfNodes());
  T1.forEdges([&] (node u, node v) {
    edgeweight w = T1.getWeight(T1.getParent(u), u);
    EXPECT_LE(lower - EPSILON, w);
    EXPECT_LE(w, upper + EPSILON);
  });
}

/* Tests the computation of resistances. */
TEST_F(RootedTreeGTest, testResistances) {
  /**
   *      3           5
   *     /           /
   *  (1/2)       (1/2)
   *   /           /
   *  0---(1/3)---2
   *   \           \
   *  (1/2)        (1)
   *     \           \
   *      3           4
   */
  RootedTree T({
    {0, 1, 1./2}, {0, 2, 1./3}, {0, 3, 1./2}, {2, 4, 1./1}, {2, 5, 1./2}
  }, 0);
  std::vector<Edge> edges = {{1, 3}, {1, 4}, {4, 5}};
  vector<edgeweight> depths_expect = {0, 2, 3, 2, 4, 5};
  vector<edgeweight> resistances_expect = {4, 6, 3};

  auto depths_actual = computeReciprocalDepths(T);
  auto resistances_actual = computeResistances(T, edges);
  EXPECT_EQ(depths_expect, depths_actual);
  EXPECT_EQ(resistances_expect, resistances_actual);
}

}
}

/** @endinternal */
#endif

