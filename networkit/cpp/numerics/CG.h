/*
 * CG.h
 *
 *  Created on: 15.06.2014
 *      Author: dhoske
 */

#ifndef CG_H_
#define CG_H_

#include "../algebraic/Matrix.h"
#include "../algebraic/Vector.h"

namespace NetworKit {

/**
 * Solves the linear system \f$Ax = b\f$ using the conjugate gradient method
 * without preconditioning and with initial value \f$(0, \dots, 0)^T\f$.
 * The returned solution \f$x\f$ fulfils \f$\frac{\Vert Ax - b\Vert}{\Vert b \Vert}
 * \leq relative\_residual\f$.
 *
 * Obviously, @a A needs to have the same number of rows as @a b and
 * @a relative_residual must be positive.
 */
Vector conjugateGradient(const Matrix& A, const Vector& b, double relative_residual);

}
#endif
