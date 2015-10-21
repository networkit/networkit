/*
 * Lamg.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: Michael
 */

#include "Lamg.h"

#include "../../properties/ConnectedComponents.h"
#include "../GaussSeidelRelaxation.h"

namespace NetworKit {


Lamg::Lamg(const double desiredResidualRed) : LinearSolver(desiredResidualRed), validSetup(false), lamgSetup(smoother) {
}

void Lamg::setup(const CSRMatrix &laplacianMatrix) {
	this->laplacianMatrix = laplacianMatrix;
	Graph G = CSRMatrix::matrixToGraph(laplacianMatrix);
	ConnectedComponents con(G);
	con.run();
	graph2Components = std::vector<index>(G.numberOfNodes());
//	components2Graph = std::vector<std::vector<index>>(con.numberOfComponents(), std::vector<index>());

	initialVectors = std::vector<Vector>(con.numberOfComponents());
	rhsVectors = std::vector<Vector>(con.numberOfComponents());

	components = std::vector<std::vector<index>>(con.numberOfComponents());
	compHierarchies = std::vector<LevelHierarchy>(con.numberOfComponents());
	compSolvers.clear();
	compStati = std::vector<LAMGSolverStatus>(con.numberOfComponents());

	// create solver for every component
	index compIdx = 0;
	for (auto component : con.getPartition().getSubsets()) {
		components[compIdx] = std::vector<index>(component.begin(), component.end());

		std::vector<std::pair<index,index>> positions;
		std::vector<double> values;
//		components2Graph[compIdx] = std::vector<index>(component.size());

		index idx = 0;
		for (node u : components[compIdx]) {
			graph2Components[u] = idx;
//			components2Graph[compIdx][idx] = u;
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

SolverStatus Lamg::solve(const Vector &rhs, Vector &result, count maxConvergenceTime, count maxIterations) {
	if (!validSetup || result.getDimension() != laplacianMatrix.numberOfColumns()
								|| rhs.getDimension() != laplacianMatrix.numberOfRows()) {
		throw std::runtime_error("No or wrong matrix is setup for given vectors.");
	}

	// solve on every component
	count maxIters = 0;
	for (index i = 0; i < components.size(); ++i) {
		for (auto element : components[i]) {
			initialVectors[i][graph2Components[element]] = result[element];
			rhsVectors[i][graph2Components[element]] = rhs[element];
		}

		compStati[i].maxIters = maxIterations;
		compStati[i].maxConvergenceTime = maxConvergenceTime;
		compSolvers[i].solve(initialVectors[i], rhsVectors[i], compStati[i]);

		for (auto element : components[i]) { // write solution back to result
			result[element] = initialVectors[i][graph2Components[element]];
		}

		maxIters = std::max(maxIters, compStati[i].numIters);
	}

	SolverStatus status;

	if (components.size() > 1) {
		status.residual = (rhs - laplacianMatrix * result).length();
		status.converged = status.residual <= tolerance;
		status.numIters = maxIters;
	} else {
		status.residual = compStati[0].residual;
		status.numIters = compStati[0].numIters;
		status.converged = compStati[0].converged;
	}

	return status;
}



} /* namespace NetworKit */
