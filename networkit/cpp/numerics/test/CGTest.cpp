/*
 * CGTest.cpp
 *
 *  Created on: 15.06.2014
 *      Author: dhoske
 */

#ifndef NOGTEST

#include <vector>
#include <string>
#include <random>
#include <functional>

#include "CGTest.h"
#include "../CG.h"
#include "../../algebraic/LaplacianMatrix.h"
#include "../../io/METISGraphReader.h"

using namespace std;

namespace NetworKit {

// Settings for CG tests
// Relative residual
static const double RESIDUAL = 1e-6;
// Graphs to test Laplace matrix on
static vector<string> GRAPHS = {"jazz", "power", "PGPgiantcompo"};
// Seed for generation of random vector
static const size_t SEED = 4342533;

// Helper function to compute relative residual
static double residual(const Matrix& A, const Vector& b, const Vector& x) {
  return (A*x - b).length() / b.length();
}

// Test on small manually generated matrix
TEST_F(CGTest, testSmall) {
  vector<pair<int, int>> idx = {{0, 0}, {0, 1}, {1, 1}};
  Matrix A(2, idx, {2, 1, 3});
  Vector b = {0, 5};
  Vector x = conjugateGradient(A, b, RESIDUAL);
  EXPECT_LE(residual(A, b, x), RESIDUAL);
}

// Test on larger Laplace matrixes (is this test too slow?)
TEST_F(CGTest, testLargeLaplace) {
  METISGraphReader reader;
  auto rand_val = bind(
        uniform_real_distribution<double>(-1.0, 1.0),
        default_random_engine(SEED)
  );

  for (const string& name : GRAPHS) {
    Graph G = reader.read("input/" + name + ".graph");
    count n = G.numberOfNodes();
    LaplacianMatrix L(G);

    // Generate random vector that sums to zero (should be in image of L)
    Vector b(G.numberOfNodes());
    double sum = 0.0;
    for (index i = 0; i < n - 1; ++i) {
      b[i] = rand_val();
      sum += b[i];
    }
    b[n - 1] = -sum;

    // Now check the CG method
    Vector x = conjugateGradient(L, b, RESIDUAL);
    EXPECT_LE(residual(L, b, x), RESIDUAL);
  }
}

}

#endif
