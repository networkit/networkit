/*
 * GaussSeidelRelaxation.h
 *
 *  Created on: 27.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef GAUSSSEIDELRELAXATION_H_
#define GAUSSSEIDELRELAXATION_H_

#include "Smoother.h"

namespace NetworKit {

/**
 * @ingroup numerics
 * Implementation of the Gauss-Seidel smoother.
 */
class GaussSeidelRelaxation : public Smoother {

private:
	double tolerance;

public:
	/**
	 * Constructs a Gauss-Seidel smoother with the given @a tolerance (default: 1e-15).
	 * @param tolerance
	 */
	GaussSeidelRelaxation(double tolerance=1e-15);

	/**
	 * Utilizes Gauss-Seidel relaxations until the given number of @a maxIterations is reached or the relative residual
	 * is below the tolerance specified in the constructor. The solver starts with @a initialGuess as intitial guess to
	 * the solution.
	 * @param A The matrix.
	 * @param b The right-hand-side.
	 * @param initialGuess
	 * @param maxIterations
	 * @return The (approximate) solution to the system.
	 */
	Vector relax(const CSRMatrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations = std::numeric_limits<count>::max()) const;

	/**
	 * Utilizes Gauss-Seidel relaxations until the given number of @a maxIterations is reached or the relative residual
	 * is below the tolerance specified in the constructor.
	 * @param A The matrix.
	 * @param b The right-hand-side.
	 * @param maxIterations
	 * @return The (approximate) solution to the system.
	 */
	Vector relax(const CSRMatrix &A, const Vector &b, const count maxIterations = std::numeric_limits<count>::max()) const;

};

} /* namespace NetworKit */

#endif /* GAUSSSEIDELRELAXATION_H_ */
