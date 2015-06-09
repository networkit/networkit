///*
// * LAMGSolver.cpp
// *
// *  Created on: 25.11.2014
// *      Author: Michael
// */
//
//#include "LAMGSolver.h"
//#include "../auxiliary/Enforce.h"
//#include <sstream>
//#include <fstream>
//#include <iostream>
//
//namespace NetworKit {
//
//LAMGSolver::LAMGSolver(double tolerance, const MultigridHierarchyBuilder<LAMGHierarchy> &builder, const Smoother &smoother) : LinearSolver(tolerance), builder(builder), smoother(smoother) {
//}
//
//void LAMGSolver::setup(const Matrix &matrix) {
//	hierarchy = builder.buildHierarchy(matrix);
//}
//
//Vector LAMGSolver::multigridCycle(const index level, const Vector &x, const Vector &bf) {
//	if (level == hierarchy.getNumLevels() - 1) { // coarsest level
////		CGStatus status;
////		status.residual = tolerance;
////		status.max_iters = 10000;
////		return solveConjugateGradient<DiagonalPreconditioner>(hierarchy.getLaplacian(level), rhs, status);
//		return smoother.relax(hierarchy.getLaplacian(level), bf, x, 1e4);
//	}
//
//	// perform pre-smoothing if necessary
//	Vector xf = smoother.relax(hierarchy.getLaplacian(level), bf, x, hierarchy.getNumPreSmooth(level));
//
//	count cycles = hierarchy.getNumMultigridCycles(level);
////	INFO("cycles = ", cycles);
//	Vector xc = hierarchy.createInitialCoarseResult(level, xf);
//	// restriction to coarser level
//	Vector bc = hierarchy.restriction(level, xf, bf);
//
//	for (count cycleCount = 0; cycleCount < cycles; ++cycleCount) {
//		// recursively solve on next coarser level with restrictedRHS
//		xc = multigridCycle(level + 1, xc, bc);
//
////		INFO(level, ": ", (hierarchy.getLaplacian(level + 1) * xc - bc).length());
//
//		// TODO: Save previous pre-relaxed iterate at this level before interpolating and relaxing it
//
//
//
//
//	}
//
//	// TODO: iterate recombination
//
//	xf = hierarchy.prolongation(level, xc, xf, bf);
//
//	// perform post-smoothing if necessary
//	xf = smoother.relax(hierarchy.getLaplacian(level), bf, xf, hierarchy.getNumPostSmooth(level));
//
//	return xf;
//}
//
////Vector LAMGSolver::cycle(const Vector &xFine, const Vector &bFine) {
////	index level = 0;
////	index k = 0;
////	float maxVisits = 0.0;
////	std::vector<count> numVisits(hierarchy.getNumLevels(), 0);
////
////	while (true) {
////		if (level == hierarchy.getNumLevels() - 1) {
////			k = level-1;
////
////		} else {
////			if (level == 0) { // finest level
////				maxVisits = 1.0;
////			} else {
////				maxVisits = hierarchy.getNumMultigridCycles(level) * numVisits[level-1];
////			}
////
////			if (numVisits[level] < maxVisits) {
////				k = level+1;
////			} else {
////				k = level-1;
////			}
////		}
////
////
////	}
////
////
////}
//
//bool LAMGSolver::solve(const Vector &rhs, Vector &result, count maxIterations) {
//	count iter = 0;
//	Matrix matrix = hierarchy.getLaplacian(0);
//	Vector residual = rhs - matrix * result;
//	Vector initial = residual;
//
//	INFO("Entering while loop to solve equation");
//
//	while (!isConverged(residual) /*&& iter < maxIterations*/) {
//		INFO("residual norm= ", residual.length());
//		result = multigridCycle(0, result, rhs);
//
//		double mean;
//		for (index i = 0; i < result.getDimension(); ++i) {
//			mean += result[i];
//		}
//		mean /= result.getDimension();
//
//		for (index i = 0; i < result.getDimension(); ++i) {
//			result[i] -= mean;
//		}
//
//		residual = rhs - matrix * result;
//		++iter;
//	}
//
//	double acf = pow(residual.length() / initial.length(), (double) 1 / iter);
//	INFO("acf=", acf);
//
//	return isConverged(residual);
//}
//
//bool LAMGSolver::solve(const std::string &graph_name, const Vector &rhs, Vector &x, count maxIterations) {
//	count iter = 0;
//	Matrix matrix = hierarchy.getLaplacian(0);
//	Vector residual = rhs - matrix * x;
//	Vector initial = residual;
//
//	INFO("Entering while loop to solve equation");
//
//	std::ofstream file("LAMG_" + graph_name + ".out");
//	Aux::enforceOpened(file);
//
//	while (!isConverged(residual) /*&& iter < maxIterations*/) {
//		INFO("residual norm= ", residual.length());
//		file << iter << " " << residual.length() << std::endl;
//		x = multigridCycle(0, x, rhs);
//
//		double mean = 0.0;
//		for (index i = 0; i < x.getDimension(); ++i) {
//			mean += x[i];
//		}
//		mean /= x.getDimension();
//
//		for (index i = 0; i < x.getDimension(); ++i) {
//			x[i] -= mean;
//		}
//
//		residual = rhs - matrix * x;
//		++iter;
//	}
//	file << iter << " " << residual.length() << std::endl;
//
//	file.close();
//
//	double acf = pow(residual.length() / initial.length(), (double) 1 / iter);
//	file.open("ACF", std::ios_base::app);
//	file << graph_name << "\t" << acf << "\n";
//	file.close();
//
//	return isConverged(residual);
//}
//
//} /* namespace NetworKit */
