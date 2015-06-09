/*
 * LaplacianPreconditioner.h
 *
 *  Created on: 04.07.2014
 *      Author: dhoske
 */

#ifndef LAPLACIANPRECONDITIONER_H_
#define LAPLACIANPRECONDITIONER_H_

#include <cstdint>
#include <utility>
#include <random>

#include "../algebraic/Matrix.h"
#include "../algebraic/Vector.h"
#include "../sdd/SDDSolver.h"

namespace NetworKit {

/**
 * @addtogroup numerics
 * @{
 */

/**
 * Preconditioner using the nearly-linear time algorithm for Laplacian matrices.
 */
template<std::uint64_t iterations>
class LaplacianPreconditioner {
public:
  LaplacianPreconditioner(const Matrix& L) : L(L) {
    assert(L.numberOfColumns() == L.numberOfRows());
    assert(SDD::isSDD(L));
  };

  /**
   * Returns the preconditioned right-hand-side \f$L^{-1}b\f$.
   */
  Vector rhs(const Vector& b) {
    SDD::SolverStatus status;
    status.max_iters = iterations;
    std::random_device rd;
    status.seed = rd();
    return SDD::solveLaplacian<SDD::StretchCycleDistribution, SDD::TrivialFlow, SDD::specialGridST>(SDD::laplacianToGraph(L), b, status);
  }

private:
  Matrix L;
};

/** @} */

}

#endif
