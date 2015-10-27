/*
 * FlowGTest.cpp
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef NOGTEST
/** @internal */

#include <memory>
#include <random>
#include <functional>

#include "../../graph/GraphGenerator.h"
#include "../Config.h"
#include "../RandomTreeGenerator.h"
#include "FlowGTest.h"

using namespace std;

namespace NetworKit {
namespace SDD {

/* Test for initial flow. */
TEST_F(FlowGTest, testInitialFlow) {
  RootedTree T({
    {0, 1, 0},
    {0, 2, 0},
    {0, 3, 0},
    {2, 4, 0},
    {2, 5, 0},
  }, 0);
  Vector b = {1, -1, -4, -1, 2, 3};
  RootedTree actual = initialFlow(T, b);
  RootedTree expected({
    {0, 1, -1},
    {0, 2, 1},
    {0, 3, -1},
    {2, 4, 2},
    {2, 5, 3},
  }, 0);
  EXPECT_EQ(expected, actual);
}

/* Basic flow test on manual graph. */
TEST_F(FlowGTest, testBasic) {
  for (int i = 1; i <= 3; ++i) {
    RootedTree T({
      {0, 1, 1},
      {0, 2, 2},
      {0, 3, 3},
      {2, 4, 4},
      {2, 5, 5},
    }, 0);

    unique_ptr<Flow> F;
    Vector b(T.numberOfNodes(), 0.0);
    if (i == 1) {
      F.reset(new TrivialFlow(T, b));
    } else if (i == 2) {
      F.reset(new LogFlow(T, b));
    } else {
      F.reset(new LCAFlow(T, b));
    }
    F->addFlow(1, 0, 10);
    EXPECT_EQ(10, F->getPotential(1, 0));
    F->addFlow(4, 0, 12);
    EXPECT_EQ(9, F->getPotential(4, 0));
    F->addFlow(5, 0, 5);
    EXPECT_EQ(9.5, F->getPotential(5, 0));
    EXPECT_EQ(11.5, F->getPotential(4, 0));
    EXPECT_EQ(8.5, F->getPotential(2, 0));
    F->addFlow(2, 0, 3);
    EXPECT_EQ(11, F->getPotential(5, 0));
    EXPECT_EQ(13, F->getPotential(4, 0));
    EXPECT_EQ(10, F->getPotential(2, 0));
    Vector actual = F->primalSolution();
    Vector expected = {0, 10, 10, 0, 13, 11};
    EXPECT_EQ(expected, actual);
  }
}

/* Compare trivial flow and logarithmic flow results on a random tree. */
TEST_F(FlowGTest, testCompareRandom) {
  static const int NODES = 1000;
  static const unsigned SEED = 12345;
  default_random_engine rand(SEED);
  auto random_node = uniform_int_distribution<node>(0, NODES - 1);
  auto random_flow = uniform_real_distribution<edgeweight>(1.0, 1000.0);

  RandomTreeGenerator gen(SEED);
  auto T = gen.makeTree(NODES, 1.0, 1000.0);
  Vector b(T.numberOfNodes(), 0.0);
  TrivialFlow F1(T, b);
  LogFlow F2(T, b);
  LCAFlow F3(T, b);

  /* Do the same random operations on both flows and compare the results */
  for (int i = 0; i < 2*NODES; ++i) {
    node n1 = random_node(rand), n2 = random_node(rand);
    double flow = random_flow(rand);
    F1.addFlow(n1, n2, flow);
    F2.addFlow(n1, n2, flow);
    F3.addFlow(n1, n2, flow);
    EXPECT_NEAR(F1.getPotential(n1, n2), F2.getPotential(n1, n2), EPSILON);
    EXPECT_NEAR(F1.getPotential(n1, n2), F3.getPotential(n1, n2), EPSILON);
  }

  auto v1 = F1.primalSolution();
  auto v2 = F2.primalSolution();
  auto v3 = F3.primalSolution();
  for (int i = 0; i < NODES; ++i) {
    EXPECT_NEAR(v1[i], v2[i], EPSILON);
    EXPECT_NEAR(v1[i], v3[i], EPSILON);
  }
}

}
}

/** @endinternal */
#endif
