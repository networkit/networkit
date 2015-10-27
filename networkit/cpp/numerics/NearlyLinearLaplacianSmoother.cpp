/*
 * NearlyLinearSDDSmoother.cpp
 *
 *  Created on: Jun 6, 2015
 *      Author: Michael
 */

#include "NearlyLinearLaplacianSmoother.h"
#include "../sdd/LaplacianSolver.h"

namespace NetworKit {

NearlyLinearLaplacianSmoother::NearlyLinearLaplacianSmoother(double tolerance) : tolerance(tolerance) {
}

Vector NearlyLinearLaplacianSmoother::relax(const CSRMatrix &A, const Vector &b, const Vector &intialGuess, count maxIterations) const {
	return relax(A, b, maxIterations);
}

Vector NearlyLinearLaplacianSmoother::relax(const CSRMatrix &A, const Vector &b, count maxIterations) const {
	Graph G = CSRMatrix::laplacianToGraph(A);
	SDD::SolverStatus status;
	status.max_iters = maxIterations;
	status.desired_residual = tolerance;

	return SDD::solveLaplacian<SDD::UniformCycleDistribution, SDD::LogFlow, SDD::minDistanceST>(G, b, status);
}

} /* namespace NetworKit */
