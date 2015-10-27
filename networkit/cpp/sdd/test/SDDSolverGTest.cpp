/*
 * SDDSolverGTest.cpp
 *
 *  Created on: May 04, 2014
 *      Author: dhoske
 */

#ifndef NOGTEST
/** @internal */

#include <atomic>
#include "SDDSolverGTest.h"

using namespace std;

namespace NetworKit {
namespace SDD {

static constexpr double TEST_RESIDUAL = 1e-9;

/* Compare two graphs (could also be moved to Graph::operator==) */
bool compare_graphs(const Graph& G1, const Graph& G2) {
  if (G1.numberOfNodes() != G2.numberOfNodes() || G1.numberOfEdges() != G2.numberOfEdges()) {
    return false;
  }

  bool result = true;
  G1.forEdges([&] (node u, node v) {
    if (G1.weight(u, v) != G2.weight(u, v)) {
      result = false;
    }
  });
  return result;
}

/* Test SDD-determination on a small manual example. */
TEST_F(SDDSolverGTest, testIsSDD) {
  vector<pair<index, index>> idx = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  Matrix m1(2, idx, {2, 1, 1, 2});
  EXPECT_TRUE(isSDD(m1));
  Matrix m2(2, idx, {2, 1, 1, 0});
  EXPECT_FALSE(isSDD(m2));
  Matrix m3(2, idx, {1, -2, -2, 1});
  EXPECT_FALSE(isSDD(m3));
}

/* Test Laplacian checker on a small manual example. */
TEST_F(SDDSolverGTest, testIsLaplacian) {
  vector<pair<index, index>> idx = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  Matrix m1(2, idx, {-1, 1, 1, -1});
  EXPECT_FALSE(isLaplacian(m1));
  Matrix m2(2, idx, {1, -1, -1, 1});
  EXPECT_TRUE(isLaplacian(m2));
}

TEST_F(SDDSolverGTest, testLaplacianToGraph) {
  vector<pair<index, index>> idx = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  Matrix m(2, idx, {2, -2, -2, 2});
  Graph G = laplacianToGraph(m);
  Graph Gexpect = {{0, 1, 2}};
  EXPECT_TRUE(compare_graphs(G, Gexpect));
}

/* Test transformation to a graph on a small manual example. */
TEST_F(SDDSolverGTest, testToGraph) {
  vector<pair<index, index>> idx = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  Vector b = {1, 1};
  Vector b_expected = {1, 1, -1, -1};
  Graph G1_expected = {{0, 1, 2}, {0, 2, 1}, {1, 3, 1}, {2, 3, 2}};
  Graph G2_expected = {{0, 2, 1}, {0, 3, 2}, {1, 2, 2}, {1, 3, 1}};

  Matrix A1(2, idx, {4, -2, -2, 4});
  Graph G1 = sddToLaplacian(A1);
  Vector b1 = sddToLaplacian(b);
  EXPECT_TRUE(compare_graphs(G1_expected, G1));
  EXPECT_EQ(b_expected,  b1);

  Matrix A2(2, idx, {4,  2, 2, 4});
  Graph G2 = sddToLaplacian(A2);
  Vector b2 = sddToLaplacian(b);
  EXPECT_TRUE(compare_graphs(G2_expected, G2));
  EXPECT_EQ(b_expected,  b2);
}

/* Test of the complete solver on a manual matrix. */
TEST_F(SDDSolverGTest, testSimple) {
  vector<pair<index, index>> idx = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  Matrix A(2, idx, {2, 1, 1, 3});
  Vector b = {0, 5};

  auto solve1 = solveSDD<UniformCycleDistribution, TrivialFlow, minDistanceST>;
  auto solve2 = solveSDD<StretchCycleDistribution, LogFlow, lowStretchST>;

  SolverStatus status;
  status.desired_residual = TEST_RESIDUAL;
  Vector x1 = solve1(A, b, status);
  EXPECT_TRUE(status.converged);
  EXPECT_EQ(uint64_t(1), status.niters);
  Vector x2 = solve2(A, b, status);
  EXPECT_TRUE(status.converged);
  EXPECT_EQ(uint64_t(0), status.niters);
  Vector x_expected = {-1.0, 2.0};
  EXPECT_EQ(x_expected, x1);
  EXPECT_ANY_THROW(solve2(Matrix(), Vector(), status));
}

}
}

/** @endinternal */
#endif

