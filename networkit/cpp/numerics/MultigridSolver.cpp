/*
 * MultigridSolver.cpp
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultigridSolver.h"
#include "../auxiliary/Enforce.h"
#include <sstream>
#include <fstream>
#include <iostream>

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
//		CGStatus status;
//		status.residual = tolerance;
//		status.max_iters = 10000;
//		return solveConjugateGradient<IdentityPreconditioner>(matrix, rhs, status);
		return smoother.relax(matrix, rhs, result, 1e4);
	}

	// presmoothing
	Vector res = smoother.relax(matrix, rhs, result, numPreSmooth);

	// restriction
	Vector restrictedVector = hierarchy.restriction(level, res, rhs);

	// coarse error correction
	Vector coarseResult(restrictedVector.getDimension());
	coarseResult = vCycle(hierarchy.getLaplacian(level + 1), level + 1, restrictedVector, coarseResult);

	// prolongation
	res = hierarchy.prolongation(level, coarseResult, res);

	// postsmooothing
	res = smoother.relax(matrix, rhs, res, numPostSmooth);

	return res;
}

Vector MultigridSolver::fCycle(const Matrix &matrix, index level, const Vector &rhs, const Vector &result) {
	if (level == hierarchy.getNumLevels() - 1) {
//		CGStatus status;
//		status.residual = tolerance;
//		status.max_iters = 10000;
//		return solveConjugateGradient<IdentityPreconditioner>(matrix, rhs, status);
		return smoother.relax(matrix, rhs, result, 1e4);
	}

	// restriction
	Vector restrictedVector = hierarchy.restriction(level, result, rhs);

	// coarse error correction (f-cycle)
	Vector coarseResult(restrictedVector.getDimension());
	coarseResult = fCycle(hierarchy.getLaplacian(level + 1), level + 1, restrictedVector, coarseResult);

	// prolongation
	Vector res = hierarchy.prolongation(level, coarseResult, result);

	// v-cycle
	res = vCycle(matrix, level, rhs, res);

	return res;
}

bool MultigridSolver::solve(const Vector &rhs, Vector &result, count maxIterations) {
	count iter = 0;
	Matrix matrix = hierarchy.getLaplacian(0);
	Vector residual = rhs - matrix * result;

	INFO("Entering while loop to solve equation");
	while (!isConverged(residual) && iter < maxIterations) {
		INFO("residual norm= ", residual.length());
		result = fCycle(matrix, 0, rhs, result);
		residual = rhs - matrix * result;

		++iter;
	}

	return isConverged(residual);
}

bool MultigridSolver::solve(const std::string &graph_name, const Vector &rhs, Vector &result, count maxIterations) {
	count iter = 0;
	Matrix matrix = hierarchy.getLaplacian(0);
	Vector residual = rhs - matrix * result;
	Vector initial = residual;

	std::ofstream file("LAMG_" + graph_name + ".out");
	Aux::enforceOpened(file);

	INFO("Entering while loop to solve equation");
	while (!isConverged(residual) /*&& iter < maxIterations*/) {
		file << iter << " " << residual.length() << std::endl;
		INFO("residual norm= ", residual.length());
		result = fCycle(matrix, 0, rhs, result);
		residual = rhs - matrix * result;

		++iter;
	}

	file << iter << " " << residual.length() << std::endl;

	file.close();

	double acf = pow(residual.length() / initial.length(), (double) 1 / iter);
	INFO("acf=", acf);

	return isConverged(residual);
}

} /* namespace NetworKit */
