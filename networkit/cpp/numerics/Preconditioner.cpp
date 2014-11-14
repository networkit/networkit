/*
 * Preconditioner.cpp
 *
 *  Created on: 05.07.2014
 *      Author: dhoske
 */

#include "Preconditioner.h"

namespace NetworKit {

Vector IdentityPreconditioner::rhs(const Vector& b) {
  return b;
}

Vector DiagonalPreconditioner::rhs(const Vector& b) {
  assert(b.getDimension() == inv_diag.getDimension());
  Vector out(b.getDimension());
  for (index i = 0; i < b.getDimension(); ++i) {
    out[i] = inv_diag[i] * b[i];
  }
  return out;
}

}
