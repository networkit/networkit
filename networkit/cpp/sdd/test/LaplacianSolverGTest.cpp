/*
 * LaplacianSolverGTest.cpp
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef NOGTEST
/** @internal */

#include <iostream>
#include <random>
#include <future>

#include "LaplacianSolverGTest.h"
#include "../../graph/GraphGenerator.h"

using namespace std;

namespace NetworKit {
namespace SDD {

namespace {
  /* Try large settings */
  const vector<string> TEST_GRAPHS = {"power.graph", "lesmis.graph"};
  const vector<string> TEST_SPARSIFY = {"tiny_01.graph"};
  const SeedType TEST_SEED = 1234;
  const double TEST_RESIDUAL = 1e-2;
}

TEST_F(LaplacianSolverGTest, testInducedSubgraph) {
  GraphGenerator gen;
  auto K6 = gen.makeCompleteGraph(6);
  auto K3 = gen.makeCompleteGraph(3);
  auto induced = inducedSubgraph(K6, {1, 3, 5});
  EXPECT_EQ(K3.numberOfNodes(), induced.numberOfNodes()) << "K_6[1,3,5] should have 3 vertices";
  EXPECT_EQ(K3.numberOfEdges(), induced.numberOfEdges()) << "K_6[1,3,5] should have 3 edges";
}


TEST_F(LaplacianSolverGTest, testGenGrid) {
  Graph Gexpected = {{0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {2, 3, 1}};
  Graph Gactual = gen2DGrid(2);
  EXPECT_EQ(Gexpected.edges(), Gactual.edges());
  Gexpected.forEdges([&] (node u, node v) {
    Gactual.hasEdge(u, v);
  });
}

TEST_F(LaplacianSolverGTest, testInducedVector) {
  Vector v = {3, -1, 4, 5, 2};
  Vector expected = {5, 3, -1};
  Vector actual = inducedVector(v, {3, 0, 1});
  EXPECT_EQ(expected, actual);
}

/* Test of complete Laplace solver test on a cycle. */
TEST_F(LaplacianSolverGTest, testSimple) {
  static constexpr int SIZE = 50;
  GraphGenerator gen;
  Graph G = gen.makeCircularGraph(SIZE);
  LaplacianMatrix L(G);

  vector<SolverFunc*> solvers = {
    solveLaplacian<UniformCycleDistribution, TrivialFlow, minDistanceST>,
    solveLaplacian<StretchCycleDistribution, LogFlow, lowStretchST>
  };
  Vector b = Vector(SIZE, 0);
  for (int i = 0; i < SIZE/2; ++i) {
    b[i] = i;
    b[i + SIZE/2] = -i;
  }

  for (auto& solver: solvers) {
    SolverStatus status;
    status.seed = TEST_SEED;
    status.desired_residual = TEST_RESIDUAL;
    status.max_iters = 100;

    auto v = solver(G, b, status);
    EXPECT_TRUE(status.converged);
    EXPECT_EQ(status.niters, 1u);
    EXPECT_NEAR(0.0, residual(L, v, b), EPSILON);
    EXPECT_NEAR(EPSILON, 2.*(SIZE - 1)/SIZE, status.stretch);

    EXPECT_ANY_THROW(solver(Graph(), Vector(), status));
    EXPECT_ANY_THROW(solver(G, Vector(), status));
  }
}

/* Tests that the flow update is correct. */
TEST_F(LaplacianSolverGTest, testFlowUpdate) {
  Graph G = {{0, 1, 1./1}, {1, 2, 1./2}, {1, 3, 1./4}, {2, 3, 1./4}};
  LaplacianMatrix L(G);
  Vector b = {1, 1, -1, -1};

  auto solver = solveLaplacian<UniformCycleDistribution, TrivialFlow, minDistanceST>;
  SolverStatus status;
  status.desired_residual = TEST_RESIDUAL;
  status.max_iters = 1;
  auto v = solver(G, b, status);
  EXPECT_TRUE(status.converged);
  EXPECT_EQ(uint64_t(1), status.niters);
  EXPECT_NEAR(0.0, residual(L, v, b), TEST_RESIDUAL);
}

/* Test of the complete solver on larger graphs (this test is used for profiling). */
TEST_F(LaplacianSolverGTest, tryLarge) {
  /* Trivial and improved solver. */
  vector<SolverFunc*> solvers = {
    solveLaplacian<StretchCycleDistribution, LCAFlow, minWeightST>,
    solveLaplacian<UniformCycleDistribution, TrivialFlow, minDistanceST>
  };

  for (index i = 0; i < solvers.size(); ++i) {
    for (const string& name : TEST_GRAPHS) {
      /* Read input */
      Graph G = readGraph("input/" + name);
      LaplacianMatrix L(G);
      Vector b = randZeroSum(G, TEST_SEED);

      /* Call solver */
      SolverStatus status;
      status.max_iters = 10000000;
      status.desired_residual = TEST_RESIDUAL;
      Vector x = solvers[i](G, b, status);
      EXPECT_TRUE(status.converged);
      EXPECT_LT(status.niters, status.max_iters);
      EXPECT_NEAR(0.0, residual(L, x, b), TEST_RESIDUAL);
    }
  }
}

}
}

/** @endinternal */
#endif
