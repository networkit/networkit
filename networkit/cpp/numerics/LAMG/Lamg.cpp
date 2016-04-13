/*
 * Lamg.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Lamg.h"

#include "../../components/ParallelConnectedComponents.h"
#include "../GaussSeidelRelaxation.h"
#include "omp.h"

namespace NetworKit {


Lamg::Lamg(const double tolerance) : LinearSolver(tolerance), validSetup(false), lamgSetup(smoother), numComponents(0) {
}

void Lamg::initializeForOneComponent() {
	compHierarchies = std::vector<LevelHierarchy>(1);
	lamgSetup.setup(laplacianMatrix, compHierarchies[0]);
	compSolvers.clear();
	compSolvers.push_back(SolverLamg(compHierarchies[0], smoother));
	validSetup = true;
}

void Lamg::setupConnected(const CSRMatrix &laplacianMatrix) {
	this->laplacianMatrix = laplacianMatrix;
	initializeForOneComponent();
	numComponents = 1;
}

void Lamg::setup(const CSRMatrix &laplacianMatrix) {
	this->laplacianMatrix = laplacianMatrix;
	Graph G = CSRMatrix::matrixToGraph(laplacianMatrix);
	ParallelConnectedComponents con(G, false);
	con.run();
	numComponents = con.numberOfComponents();
	if (numComponents == 1) {
		initializeForOneComponent();
	} else {
		graph2Components = std::vector<index>(G.numberOfNodes());

		initialVectors = std::vector<Vector>(numComponents);
		rhsVectors = std::vector<Vector>(numComponents);

		components = std::vector<std::vector<index>>(numComponents);
		compHierarchies = std::vector<LevelHierarchy>(numComponents);
		compSolvers.clear();
		compStati = std::vector<LAMGSolverStatus>(numComponents);

		// create solver for every component
		index compIdx = 0;
		for (auto component : con.getPartition().getSubsets()) {
			components[compIdx] = std::vector<index>(component.begin(), component.end());

			std::vector<std::pair<index,index>> positions;
			std::vector<double> values;

			index idx = 0;
			for (node u : components[compIdx]) {
				graph2Components[u] = idx;
				idx++;
			}

			for (node u : components[compIdx]) {
				G.forNeighborsOf(u, [&](node v, edgeweight w) {
					positions.push_back(std::make_pair(graph2Components[u], graph2Components[v]));
					values.push_back(w);
				});
			}
			CSRMatrix compMatrix(component.size(), component.size(), positions, values);
			initialVectors[compIdx] = Vector(component.size());
			rhsVectors[compIdx] = Vector(component.size());
			lamgSetup.setup(compMatrix, compHierarchies[compIdx]);
			compSolvers.push_back(SolverLamg(compHierarchies[compIdx], smoother));
			LAMGSolverStatus status;
			status.desiredResidualReduction = tolerance * component.size() / G.numberOfNodes();
			compStati[compIdx] = status;

			compIdx++;
		}

		validSetup = true;
	}
}

SolverStatus Lamg::solve(const Vector &rhs, Vector &result, count maxConvergenceTime, count maxIterations) {
	if (!validSetup || result.getDimension() != laplacianMatrix.numberOfColumns()
			|| rhs.getDimension() != laplacianMatrix.numberOfRows()) {
		throw std::runtime_error("No or wrong matrix is setup for given vectors.");
	}

	SolverStatus status;

	if (numComponents == 1) {
		LAMGSolverStatus stat;
		stat.desiredResidualReduction = tolerance * rhs.length() / (laplacianMatrix * result - rhs).length();
		stat.maxIters = maxIterations;
		stat.maxConvergenceTime = maxConvergenceTime;
		compSolvers[0].solve(result, rhs, stat);

		status.residual = stat.residual;
		status.numIters = stat.numIters;
		status.converged = stat.converged;
	} else {
		// solve on every component
		count maxIters = 0;
		for (index i = 0; i < components.size(); ++i) {
			for (auto element : components[i]) {
				initialVectors[i][graph2Components[element]] = result[element];
				rhsVectors[i][graph2Components[element]] = rhs[element];
			}

			double resReduction = tolerance * rhsVectors[i].length() / (compHierarchies[i].at(0).getLaplacian() * initialVectors[i] - rhsVectors[i]).length();
			compStati[i].desiredResidualReduction = resReduction * components[i].size() / laplacianMatrix.numberOfRows();
			compStati[i].maxIters = maxIterations;
			compStati[i].maxConvergenceTime = maxConvergenceTime;
			compSolvers[i].solve(initialVectors[i], rhsVectors[i], compStati[i]);

			for (auto element : components[i]) { // write solution back to result
				result[element] = initialVectors[i][graph2Components[element]];
			}

			maxIters = std::max(maxIters, compStati[i].numIters);
		}

		status.residual = (rhs - laplacianMatrix * result).length();
		status.converged = status.residual <= tolerance;
		status.numIters = maxIters;
	}

	return status;
}

void Lamg::parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results, count maxConvergenceTime, count maxIterations) {
	if (numComponents == 1) {
		assert(rhs.size() == results.size());
		const index numThreads = omp_get_max_threads();
		if (compSolvers.size() != numThreads) {
			compSolvers.clear();

			for (index i = 0; i < (index) numThreads; ++i) {
				compSolvers.push_back(SolverLamg(compHierarchies[0], smoother));
			}
		}

		bool nested = omp_get_nested();
		if (nested) omp_set_nested(false);

#pragma omp parallel for
		for (index i = 0; i < rhs.size(); ++i) {
			index threadId = omp_get_thread_num();
			LAMGSolverStatus stat;
			stat.desiredResidualReduction = tolerance * rhs[i].length() / (laplacianMatrix * results[i] - rhs[i]).length();
			stat.maxIters = maxIterations;
			stat.maxConvergenceTime = maxConvergenceTime;
			compSolvers[threadId].solve(results[i], rhs[i], stat);
		}

		if (nested) omp_set_nested(true);
	}
}



} /* namespace NetworKit */
