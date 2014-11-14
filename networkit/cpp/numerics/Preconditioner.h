/*
 * Preconditioner.h
 *
 *  Created on: 04.07.2014
 *      Author: dhoske
 */

#ifndef PRECONDITIONER_H_
#define PRECONDITIONER_H_

#include <cstdint>
#include <utility>
#include <random>

#include "../algebraic/Matrix.h"
#include "../algebraic/Vector.h"


namespace NetworKit {

/**
 * @addtogroup numerics
 * @{
 */

/**
 * Simple preconditioner that returns the given vector unchanged.
 */
class IdentityPreconditioner {
public:
  /**
   * Constructs an identity preconditioner for the matrix A.
   */
  template<typename Matrix>
  IdentityPreconditioner(const Matrix& A) {
  }

  /**
   * Returns the preconditioned right-hand-side \f$P(b) = b\f$.
   */
  Vector rhs(const Vector& b);
};

/**
 * Simple preconditioner that approximates the matrix by a
 * diagonal matrix.
 */
class DiagonalPreconditioner {
public:
  /**
   * Constructs a diagonal preconditioner for the matrix A.
   */
  template<typename Matrix>
  DiagonalPreconditioner(const Matrix& A)
      : inv_diag(A.numberOfRows()) {
    assert(A.numberOfColumns() == A.numberOfRows());

    // Diagonal preconditioner just needs to store the inverse diagonal of A
    for (index i = 0; i < A.numberOfRows(); ++i) {
      edgeweight val = A(i, i);
      if (val) {
        inv_diag[i] = 1.0 / A(i, i);
      } else {
        inv_diag[i] = 0.0;
      }
    }
  };

  /**
   * Returns the preconditioned right-hand-side \f$P(b) = D(A)^{-1}b\f$.
   */
  Vector rhs(const Vector& b);

private:
  Vector inv_diag;
};

/** @} */

}

#endif
