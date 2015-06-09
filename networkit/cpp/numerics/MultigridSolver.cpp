///*
// * MultigridSolver.cpp
// *
// *  Created on: 31.10.2014
// *      Author: Michael Wegner (michael.wegner@student.kit.edu)
// */
//
//#include "MultigridSolver.h"
//#include "../auxiliary/Enforce.h"
//#include <sstream>
//#include <fstream>
//#include <iostream>
//
//namespace NetworKit {
//
//MultigridSolver::MultigridSolver(const double tolerance, const count numPreSmooth, const count numPostSmooth,
//		const LAMG &builder, const Smoother &smoother) : LinearSolver(tolerance), numPreSmooth(numPreSmooth), numPostSmooth(numPostSmooth), builder(builder), smoother(smoother) {
//}
//
//void MultigridSolver::setup(const Matrix &matrix) {
//	hierarchy = builder.buildHierarchy(matrix);
//}
//
//Vector MultigridSolver::restrictResidual(const Matrix &matrix, const Matrix &interpolationMatrix, const Vector &result, const Vector &rhs) {
//	Vector residual = rhs - matrix * result;
//	return Matrix::mTvMultiply(interpolationMatrix, residual);
//}
//
//void MultigridSolver::interpolateError(const Matrix &interpolationMatrix, const Vector &coarse_error, Vector &result) {
//	result += interpolationMatrix * coarse_error;
//}
//
//Vector MultigridSolver::vCycle(const Matrix &matrix, index level, const Vector &bFine, const Vector &xFine) {
//	if (level == hierarchy.getNumLevels() - 1) {
////		CGStatus status;
////		status.residual = tolerance;
////		status.max_iters = 10000;
////		return solveConjugateGradient<IdentityPreconditioner>(matrix, rhs, status);
//		return smoother.relax(matrix, bFine, xFine, 1e4);
//	}
//
//	// presmoothing
//	Vector newXFine = xFine;
//	if (hierarchy.getNumPreSmooth(level) > 0) {
//		newXFine = smoother.relax(matrix, bFine, xFine, numPreSmooth);
//	}
//
//	// restriction
//	Vector bCoarse = hierarchy.restriction(level, newXFine, bFine);
//
//	// coarse error correction
//	Vector xCoarse = hierarchy.createInitialCoarseResult(level, xFine);
//	xCoarse = vCycle(hierarchy.getLaplacian(level + 1), level + 1, bCoarse, xCoarse);
//
//	// prolongation
//	newXFine = hierarchy.prolongation(level, xCoarse, newXFine, bFine);
//
//	// postsmooothing
//	if (hierarchy.getNumPostSmooth(level) > 0) {
//		newXFine = smoother.relax(matrix, bFine, newXFine, numPostSmooth);
//	}
//
//	return newXFine;
//}
//
//Vector MultigridSolver::fCycle(const Matrix &matrix, index level, const Vector &bFine, const Vector &xFine) {
//	if (level == hierarchy.getNumLevels() - 1) {
////		CGStatus status;
////		status.residual = tolerance;
////		status.max_iters = 10000;
////		return solveConjugateGradient<IdentityPreconditioner>(matrix, rhs, status);
//		return smoother.relax(matrix, bFine, xFine, 1e4);
//	}
//
//	// restriction
//	Vector bCoarse = hierarchy.restriction(level, xFine, bFine);
//
//	// coarse error correction (f-cycle)
//	Vector xCoarse = hierarchy.createInitialCoarseResult(level, xFine);
//	xCoarse = fCycle(hierarchy.getLaplacian(level + 1), level + 1, bCoarse, xCoarse);
//
//	// prolongation
//	Vector newXFine = hierarchy.prolongation(level, xCoarse, xFine, bFine);
//
////	INFO("level= ", level, " r = ", (bFine - matrix * newXFine).length());
//
//	// v-cycle
//	newXFine = vCycle(matrix, level, bFine, newXFine);
//
////	INFO("level= ", level, " r = ", (bFine - matrix * newXFine).length());
//
//	return newXFine;
//}
//
//bool MultigridSolver::solve(const Vector &rhs, Vector &result, count maxIterations) {
//	count iter = 0;
//	Matrix matrix = hierarchy.getLaplacian(0);
//	Vector residual = rhs - matrix * result;
//
//	INFO("Entering while loop to solve equation");
//	while (!isConverged(residual) && iter < maxIterations) {
//		INFO("residual norm= ", residual.length());
//		result = fCycle(matrix, 0, rhs, result);
//		residual = rhs - matrix * result;
//
//		++iter;
//	}
//
//	return isConverged(residual);
//}
//
//bool MultigridSolver::solve(const std::string &graph_name, const Vector &rhs, Vector &result, count maxIterations) {
//	count iter = 0;
//	Matrix matrix = hierarchy.getLaplacian(0);
//	Vector residual = rhs - matrix * result;
//	Vector initial = residual;
//
//	std::ofstream file("LAMG_" + graph_name + ".out");
//	Aux::enforceOpened(file);
//
//	INFO("Entering while loop to solve equation");
//	while (initial.length() / residual.length() < 1e8 /*&& iter < maxIterations*/) {
//		file << iter << " " << residual.length() << std::endl;
//		INFO("residual norm= ", residual.length());
//
//		if (hierarchy.getType(0) == LAMGHierarchy::ELIMINATION) {
//			if (hierarchy.getLaplacian(1).numberOfRows() == 1) {
//				result = Vector(result.getDimension(), 0);
//				break;
//			} else {
//				Vector bCoarse = hierarchy.restriction(0, result, rhs);
//				Vector xCoarse = hierarchy.createInitialCoarseResult(0, result);
//
//				xCoarse = fCycle(hierarchy.getLaplacian(1), 1, bCoarse, xCoarse);
//
//				double mean = xCoarse.mean();
//				xCoarse -= mean;
//
//				result = hierarchy.prolongation(0, xCoarse, result, rhs);
//			}
//		} else {
//			result = fCycle(matrix, 0, rhs, result);
//			double mean = result.mean();
//			result -= mean;
//		}
//
//
//		//result = fCycle(matrix, 0, rhs, result);
//		residual = rhs - matrix * result;
//
//		++iter;
//	}
//
//	file << iter << " " << residual.length() << std::endl;
//
//	file.close();
//
//	double acf = pow(residual.length() / initial.length(), (double) 1 / iter);
//
//	file.open("ACF", std::ios_base::app);
//	file << graph_name << "\t" << acf << "\n";
//	file.close();
//
//	return isConverged(residual);
//}
//
//
//} /* namespace NetworKit */
