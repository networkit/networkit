/*
 * LAMG.cpp
 *
 *  Created on: 12.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LAMG.h"

namespace NetworKit {

std::vector<bool> LAMG::lowDegreeNodes(const Matrix &matrix) {
	std::vector<bool> eliminate(matrix.numberOfRows(), true);
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (matrix.nnzInRow(i) <= 4 && eliminate[i]) {
			matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight w){
				eliminate[j] = false;
			});
		}
	}

	return eliminate;
}

std::vector<Vector> LAMG::computeTestVectors(const Matrix &matrix, const count numberOfVectors) const {
	std::vector<Vector> testVectors(numberOfVectors, Vector(matrix.numberOfColumns()));
	Vector zeroVector(matrix.numberOfColumns(), 0.0);
	for (count i = 0; i < numberOfVectors; ++i) {
		for (count j = 0; j < matrix.numberOfColumns(); ++j) {
			testVectors[i][j] = Aux::Random::real(-1.0, 1.0);
		}

		// do 3 GS sweeps on the system Av = 0
		testVectors[i] = smoother.relax(matrix, zeroVector, testVectors[i], 3);
	}

	return testVectors;
}

Matrix LAMG::computeAffinityMatrix(const Matrix &matrix, const std::vector<Vector> &testVectors) const {
	assert(testVectors.size() > 0);
	Matrix C(matrix.numberOfRows(), matrix.numberOfColumns());
	matrix.forNonZeroElementsInRowOrder([&](index i, index j, edgeweight value){
		if (j < i) continue; // symmetric
		double ij = 0.0;
		double ii = 0.0;
		double jj = 0.0;
		for (count k = 0; k < testVectors.size(); ++k) {
			ii += testVectors[k][i] * testVectors[k][i];
			jj += testVectors[k][j] * testVectors[k][j];
			ij += testVectors[k][i] * testVectors[k][j];
		}

		double weight = (ij * ij) / (ii*ii * jj*jj);

		C.setValue(i, j, weight);
		C.setValue(j, i, weight); // symmetric
	});

	return C;
}

bool LAMG::addEliminationLevel(Matrix &matrix, LAMGHierarchy &hierarchy) {
	bool nodesEliminated = false;

	Matrix P = DiagonalMatrix(matrix.numberOfRows());
	while (matrix.numberOfRows() > 1) {
		std::vector<bool> eliminate = lowDegreeNodes(matrix);
		if (eliminate.size() < 0.01 * matrix.numberOfRows()) {
			if (nodesEliminated) { // add level if nodes have been eliminated
				hierarchy.addLevel(matrix, P, LAMGHierarchy::ELIMINATION);
			}
			return nodesEliminated;
		}

		count numFNodes = eliminate.size();
		count numCNodes = matrix.numberOfRows() - numFNodes;

		// calculate matrix S = {(-A_FC^T * A_FF^-1)^T, Ic} and permutation matrix M
		Matrix S(matrix.numberOfRows(), numCNodes);
		Matrix M(matrix.numberOfRows());

		index fIndex = 0;
		index cIndex = numFNodes;
		for (index i = 0; i < matrix.numberOfRows(); ++i) {
			if (eliminate[i]) {
				matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight value){ // (-A_FC^T * A_FF^-1)^T in row i
					S.setValue(fIndex, j, -value / matrix.nnzInRow(i));
				});

				M.setValue(i, fIndex, 1);
				fIndex++;
			} else {
				S.setValue(cIndex, cIndex - numFNodes, 1);
				M.setValue(i, cIndex, 1);
				cIndex++;
			}
		}

		// compute P
		Matrix P_iteration = M * S;

		// compute coarsened matrix and add level
		matrix = Matrix::mTmMultiply(P_iteration, matrix * P_iteration);
		P = P * P_iteration;
	}

	if (nodesEliminated) { // add level if nodes have been eliminated
		hierarchy.addLevel(matrix, P);
	}

	return nodesEliminated;
}

bool LAMG::addAggregationLevel(Matrix &matrix, LAMGHierarchy &hierarchy) {
	// TODO
	return false;
}

LAMGHierarchy LAMG::buildHierarchy(const Graph &graph) {
	return LAMGHierarchy();
}

LAMGHierarchy LAMG::buildHierarchy(const Matrix &matrix) {
	return LAMGHierarchy();
}

} /* namespace NetworKit */
