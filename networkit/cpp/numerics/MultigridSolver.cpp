/*
 * MultigridSolver.cpp
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultigridSolver.h"

namespace NetworKit {


MultigridSolver::MultigridSolver(const double tolerance, const count numPreSmooth, const count numPostSmooth,
		const MultigridHierarchyBuilder &builder, const Smoother &smoother) : LinearSolver(tolerance), numPreSmooth(numPreSmooth), numPostSmooth(numPostSmooth), builder(builder), smoother(smoother) {
}

void MultigridSolver::setup(const Matrix &matrix) {
	hierarchy = builder.buildHierarchy(matrix);
}

Vector MultigridSolver::restrictResidual(const Matrix &matrix, const Matrix &interpolationMatrix, const Vector &result, const Vector &rhs) {
	Vector residual = rhs - matrix * result;
	return Matrix::mTvMultiply(interpolationMatrix, residual);
}

void MultigridSolver::interpolateError(const Matrix &interpolationMatrix, const Vector &coarse_error, Vector &result) {
	result += interpolationMatrix * coarse_error;
}

Vector MultigridSolver::vCycle(const Matrix &matrix, index level, const Vector &rhs, const Vector &result) {
	if (level == hierarchy.getNumLevels() - 1) {
		// TODO: direct solution
		CGStatus status;
		status.residual = tolerance;
		status.max_iters = 10000;
		return solveConjugateGradient<IdentityPreconditioner>(matrix, rhs, status);
	}

	// presmoothing
	Vector res = smoother.relax(matrix, rhs, result, numPreSmooth);

	// restriction
	Matrix interpolationMatrix = hierarchy.getInterpolationMatrix(level);
	Vector restricted_residual = restrictResidual(matrix, interpolationMatrix, res, rhs);

	// coarse error correction
	Vector coarse_error(interpolationMatrix.numberOfColumns());
	coarse_error = vCycle(hierarchy.getLaplacian(level + 1), level + 1, restricted_residual, coarse_error);

	// prolongation
	interpolateError(interpolationMatrix, coarse_error, res);

	// postsmooothing
	res = smoother.relax(matrix, rhs, res, numPostSmooth);

	return res;
}

Vector MultigridSolver::fCycle(const Matrix &matrix, index level, const Vector &rhs, const Vector &result) {
	TRACE("entering fCycle at level ", level);
	if (level == hierarchy.getNumLevels() - 1) {
		// TODO: direct solution
		CGStatus status;
		status.residual = tolerance;
		status.max_iters = 10000;
		return solveConjugateGradient<IdentityPreconditioner>(matrix, rhs, status);
	}

	// presmoothing
	TRACE("smoothing on level ", level);
	Vector res = smoother.relax(matrix, rhs, result, numPreSmooth);

	// restriction
	TRACE("restriction on level", level);
	Matrix interpolationMatrix = hierarchy.getInterpolationMatrix(level);
	Vector restricted_residual = restrictResidual(matrix, interpolationMatrix, res, rhs);

	// coarse error correction (f-cycle)
	TRACE("coarse error correction from ", level, " to ", level +1);
	Vector coarse_error(interpolationMatrix.numberOfColumns());
	coarse_error = fCycle(hierarchy.getLaplacian(level + 1), level + 1, restricted_residual, coarse_error);

	// prolongation
	TRACE("prolongation on level", level);
	interpolateError(interpolationMatrix, coarse_error, res);

	// v-cycle
	res = vCycle(matrix, level, rhs, res);

	return result;
}

bool MultigridSolver::solve(const Vector &rhs, Vector &result, count maxIterations) {
	count iter = 0;
	Matrix matrix = hierarchy.getLaplacian(0);
	Vector residual = rhs - matrix * result;

	while (!isConverged(residual) && iter < maxIterations) {
		result = fCycle(matrix, 0, rhs, result);
		residual = rhs - matrix * result;
		++iter;
	}

	return isConverged(residual);
}

} /* namespace NetworKit */
