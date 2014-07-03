/*
 * CG.cpp
 *
 *  Created on: 15.06.2014
 *      Author: dhoske
 */

#include "CG.h"

namespace NetworKit {

Vector conjugateGradient(const Matrix& A, const Vector& b, double relative_residual) {
  assert(A.numberOfRows() == b.getDimension());
  assert(relative_residual > 0);

  // Absolute residual to achieve
  double residual = relative_residual * b.length();

  // Main loop. See: http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
  Vector x(A.numberOfColumns(), 0.0);
  Vector rest = b - A*x;
  Vector conjugate_dir = rest;
  double sqr_rest = rest.lengthSqr();
  while (rest.length() > residual) {
    Vector rest_dir = A * conjugate_dir;

    // dot() or innerProduct() does not exist???
    double step = sqr_rest / (rest_dir.transpose() * conjugate_dir);
    x    += step * conjugate_dir;
    rest -= step * rest_dir;

    double new_sqr_rest = rest.lengthSqr();
    conjugate_dir = new_sqr_rest / sqr_rest * conjugate_dir + rest;
    sqr_rest = new_sqr_rest;
  }

  return x;
}

}
