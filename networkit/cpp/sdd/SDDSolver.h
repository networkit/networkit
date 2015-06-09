/*
 * SDDSolver.h
 *
 *  Created on: May 04, 2014
 *      Author: dhoske
 */

#ifndef SDD_SOLVER_TREE_H_
#define SDD_SOLVER_TREE_H_

#include "../graph/Graph.h"
#include "../algebraic/Matrix.h"
#include "Config.h"
#include "LaplacianSolver.h"

namespace NetworKit {

namespace SDD {

/** @addtogroup sdd
 *  @{ */

/**
 * Tests whether the matrix @a A is symmetric and diagonal dominant.
 */
bool isSDD(const Matrix& A);


/**
 * Tests whether the matrix @a A is a Laplacian matrix, i.e. it is symmetric,
 * the off-diagonal-entries are non-positive and the sum of entries in each row is zero.
 */
bool isLaplacian(const Matrix& A);

/**
 * Converts a Laplacian matrix to the corresponding graph.
 */
Graph laplacianToGraph(const Matrix& A);

/**
 * Computes the graph @a G whose Laplacian corresponds to @a A when
 * a SDD-system is transformed into a Laplacian system
 */
Graph sddToLaplacian(const Matrix& A);

/**
 * Computes the vector @a x corresponding to @a b when
 * a SDD-system is transformed into a Laplacian system.
 */
Vector sddToLaplacian(const Vector& b);

/**
 * Solves the SDD-system \f$Ax = b\f$.
 *
 * The options in @a status are used both as input of the solver and
 * as output for its status information. In particular, you can set
 * the desired residual and the maximum number of iterations in @a status.
 *
 * @tparam TCycle distribution on basis cycle to use
 * @tparam TFlow data structure to use to store the flow
 * on the spanning tree edges
 * @tparam STAlgo algorithm to compute the spanning tree
 * that is used as cycle basis
 */
template<typename TCycle = UniformCycleDistribution, typename TFlow = TrivialFlow,
         SpanningTreeAlgo STAlgo = minDistanceST>
Vector solveSDD(const Matrix& A, const Vector& b, SolverStatus& status) {
  if (A.numberOfRows() == 0) {
    throw std::invalid_argument("empty matrix disallowed");
  }
  if (!isSDD(A)) {
    throw std::invalid_argument("Expected SDD matrix");
  }
  if (A.numberOfColumns() != b.getDimension()) {
    throw std::invalid_argument("b should have the same size as A");
  }
  assert(status.desired_residual >= 0.);
  count n = b.getDimension();

  /* Just pass the problem to the Laplace solver */
  Graph G = sddToLaplacian(A);
  Vector y = sddToLaplacian(b);
  Vector x = solveLaplacian<TCycle, TFlow, STAlgo>(G, y, status);

  /* Translate solution back */
  Vector xorig(x.getDimension() / 2);
  xorig.parallelForElements([&] (index i, double value) {
    xorig[i] = (x[i] - x[i + n]) / 2.0;
  });
  return xorig;
}

/** @} */

/**
 * Tests whether the matrix @a A is symmetric. Note that all
 * of our matrices are symmetric at the moment
 */
inline bool isSymmetric(const Matrix& A) {
  bool output = true;
  A.forNonZeroElementsInRowOrder([&] (index i, index j, edgeweight w) {
    if (A(j, i) != w) {
      output = false;
    }
  });
  return output;
}

}
}

#endif
