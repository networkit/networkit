/*
 * ConjugateGradient.h
 *
 *  Created on: 15.06.2014
 *      Author: Daniel Hoske and Michael Wegner
 */

#ifndef CONJUGATE_GRADIENT_H_
#define CONJUGATE_GRADIENT_H_

#include <cstdint>
#include <utility>

#include "LinearSolver.h"
#include "../algebraic/Vector.h"
#include "../algebraic/CSRMatrix.h"

namespace NetworKit {

/**
 * @ingroup numerics
 * Implementation of Conjugate Gradient.
 */
template<class Preconditioner>
class ConjugateGradient : public LinearSolver {
public:
	ConjugateGradient(double tolerance = 1e-5) : LinearSolver(tolerance), matrix(CSRMatrix()) {}

	void setup(const CSRMatrix& matrix) {
		this->matrix = matrix;
		precond = Preconditioner(matrix);
	}

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
	SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000, count maxIterations = std::numeric_limits<count>::max());

	/**
	 * Solves the linear systems in parallel.
	 * @param rhs
	 * @param results
	 * @param maxConvergenceTime
	 * @param maxIterations
	 */
	void parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results, count maxConvergenceTime = 5 * 60 * 1000, count maxIterations = std::numeric_limits<count>::max());

private:
	CSRMatrix matrix;
	Preconditioner precond;

};

template<class Preconditioner>
SolverStatus ConjugateGradient<Preconditioner>::solve(const Vector& rhs, Vector& result, count maxConvergenceTime, count maxIterations) {
	assert(matrix.numberOfRows() == rhs.getDimension());

	// Absolute residual to achieve
	double sqr_desired_residual = tolerance * tolerance * (rhs.length() * rhs.length());

	// Main loop. See: http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
	Vector residual_dir = rhs - matrix*result;
	Vector conjugate_dir = precond.rhs(residual_dir);
	double sqr_residual = Vector::innerProduct(residual_dir, residual_dir);
	double sqr_residual_precond = Vector::innerProduct(residual_dir, conjugate_dir);

	count niters = 0;
	Vector tmp, residual_precond;
	while (sqr_residual > sqr_desired_residual) {
		niters++;
		if (niters > maxIterations) {
			break;
		}

		tmp = matrix * conjugate_dir;
		double step = sqr_residual_precond / Vector::innerProduct(conjugate_dir, tmp);
		result += step * conjugate_dir;
		residual_dir -= step * tmp;
		sqr_residual = Vector::innerProduct(residual_dir, residual_dir);

		residual_precond = precond.rhs(residual_dir);
		double new_sqr_residual_precond = Vector::innerProduct(residual_dir, residual_precond);
		conjugate_dir = (new_sqr_residual_precond / sqr_residual_precond) * conjugate_dir + residual_precond;
		sqr_residual_precond = new_sqr_residual_precond;
	}

	SolverStatus status;
	status.numIters = niters;
	status.residual = (rhs - matrix*result).length();
	status.converged = status.residual / rhs.length() <= tolerance;

	return status;
}

template<class Preconditioner>
void ConjugateGradient<Preconditioner>::parallelSolve(const std::vector<Vector>& rhs, std::vector<Vector>& results, count maxConvergenceTime, count maxIterations) {
#pragma omp parallel for
	for (index i = 0; i < rhs.size(); ++i) {
		this->solve(rhs[i], results[i], maxConvergenceTime, maxIterations);
	}
}






} /* namespace NetworKit */

#endif /* CONJUGATE_GRADIENT_H_ */
