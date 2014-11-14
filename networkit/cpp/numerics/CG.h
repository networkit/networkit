/*
 * CG.h
 *
 *  Created on: 15.06.2014
 *      Author: dhoske
 */

#ifndef CG_H_
#define CG_H_

#include <cstdint>
#include <utility>

#include "../algebraic/Matrix.h"
#include "../algebraic/Vector.h"

namespace NetworKit {

/**
 * @defgroup numerics Numerics
 *
 * Various numerical algorithms.
 * @{
 */

/**
 * Configuration options and status of the conjugate gradient algorithm.
 */
struct CGStatus {
  /** Input: Desired relative residual. */
  double residual = 0.0;
  /** Input: Maximum number of iterations. */
  std::uint64_t max_iters = std::numeric_limits<uint64_t>::max();

  /** Output: Has algorithm converged? */
  bool converged = false;

  /** Input: Enable tracing the residuals. */
  bool debug_enable_residuals = false;
  /** Output: Trace of the residuals. */
  std::vector<double> debug_residuals;
};

/**
 * Solves the linear system \f$Ax = b\f$ using the conjugate gradient method
 * with a given preconditioner and with initial value \f$(0, \dots, 0)^T\f$.
 * We the return the solution \f$x\f$. The solution \f$x\f$ fulfils
 * \f$\frac{\Vert Ax - b\Vert}{\Vert b \Vert} \leq relative\_residual\f$ if the
 * algorithm has converged.
 *
 * Obviously, @a A needs to have the same number of rows as @a b and
 * @a status.residual must be nonnegative. You may also request that the algorithm
 * does not run for more than @a status.max_iters iterations.
 */
template<typename Preconditioner, typename Matrix>
Vector solveConjugateGradient(const Matrix& A, const Vector& b, CGStatus& status) {
  assert(A.numberOfRows() == b.getDimension());
  assert(status.residual > 0);
  Preconditioner precond(A);

  // Absolute residual to achieve
  double sqr_desired_residual = status.residual * status.residual * (b.length() * b.length());

  // Main loop. See: http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
  Vector x(A.numberOfColumns(), 0.0);
  Vector residual_dir = b - A*x;
  Vector conjugate_dir = precond.rhs(residual_dir);
  double sqr_residual = residual_dir.transpose() * residual_dir;
  double sqr_residual_precond = residual_dir.transpose() * conjugate_dir;

  count niters = 0;
  Vector tmp, residual_precond;
  while (sqr_residual > sqr_desired_residual) {
    if (status.debug_enable_residuals) {
      status.debug_residuals.emplace_back(sqrt(sqr_residual)/b.length());
    }

    niters++;
    if (niters > status.max_iters) {
      return x;
    }

    // dot() or innerProduct() did not exist before???
    tmp = A * conjugate_dir;
    double step = sqr_residual_precond / (tmp.transpose() * conjugate_dir);
    x += step * conjugate_dir;
    residual_dir -= step * tmp;
    sqr_residual = residual_dir.length() * residual_dir.length();

    residual_precond = precond.rhs(residual_dir);
    double new_sqr_residual_precond = residual_dir.transpose() * residual_precond;
    conjugate_dir = (new_sqr_residual_precond / sqr_residual_precond) * conjugate_dir + residual_precond;
    sqr_residual_precond = new_sqr_residual_precond;
  }

  status.converged = true;
  return x;
}

/** @} */

}
#endif
