/*
 * LAMG.cpp
 *
 *  Created on: 12.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LAMG.h"

namespace NetworKit {

LAMG::LAMG(double guard, double cycleIndex) : guard(guard), cycleIndex(cycleIndex) {
}

std::vector<bool> LAMG::lowDegreeNodes(const Matrix &matrix) const {
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

std::vector<Vector> LAMG::computeTestVectors(const Matrix &matrix, const count numOfVectors) const {
	std::vector<Vector> testVectors(numOfVectors, Vector(matrix.numberOfColumns()));
	Vector zeroVector(matrix.numberOfColumns(), 0.0);
	for (count i = 0; i < numOfVectors; ++i) {
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
		if (j >= i) {  // symmetric
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
		}
	});

	return C;
}

void LAMG::addHighDegreeSeedNodes(const Matrix &matrix, std::vector<int64_t> &status) const {
	std::vector<count> degree(matrix.numberOfRows(), 0);
#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		degree[i] = matrix.nnzInRow(i);
	}

	std::sort(degree.begin(), degree.end());
	count median = degree[floor(matrix.numberOfRows() / 2)];

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (matrix.nnzInRow(i) >= 8 * median) {
			status[i] = 0;
		}
	}
}

std::vector<std::vector<index>> LAMG::computeStrongNeighbors(const Matrix &affinityMatrix, const double delta) const {
	std::vector<std::vector<index>> strongNeighbors(affinityMatrix.numberOfRows());
	std::vector<double> maxNeighbor(affinityMatrix.numberOfRows(), std::numeric_limits<double>::min());

	affinityMatrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) { // determine the highest affinity neighbor of each node
		if (value > maxNeighbor[i]) {
			maxNeighbor[i] = value;
		}
	});


	affinityMatrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) {
		if (value >= delta * std::max(maxNeighbor[i], maxNeighbor[j])) {
			strongNeighbors[i].push_back(j);
		}
	});

	return strongNeighbors;
}

int64_t LAMG::findBestSeed(const Matrix &affinityMatrix, const std::vector<index> &strongNeighbors, const std::vector<int64_t> &status, const index u) const {
	double maxAffinity = std::numeric_limits<double>::min();
	int64_t s = -1;
	for (index i = 0; i < strongNeighbors.size(); ++i) {
		index v = strongNeighbors[i];
		if (status[v] <= 0) { // neighbors is seed or undecided
			if (affinityMatrix(u, v) > maxAffinity) {
				s = v;
				maxAffinity = affinityMatrix(u, v);
			}
		}
	}

	return s;
}

void LAMG::aggregationStage(const Matrix &matrix, count &numCoarseNodes, const Matrix &affinityMatrix, std::vector<Vector> &testVectors, std::vector<int64_t> &status, std::vector<count> &aggregateSize, double delta) const {
	std::vector<std::vector<index>> strongNeighbors = computeStrongNeighbors(affinityMatrix, delta);
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (status[i] == -1 && strongNeighbors[i].size() > 0) { // undecided nodes with strong neighbors
			int64_t s = findBestSeed(affinityMatrix, strongNeighbors[i], status, i);
			if (s != -1) {
				status[s] = 0; // s becomes seed
				status[i] = s; // i's seed is s
				numCoarseNodes--;

				for (index k = 0; k < testVectors.size(); ++k) { // update test vectors
					testVectors[k][i] = testVectors[k][s];
				}

				aggregateSize[s]++;
				aggregateSize[i] = aggregateSize[s];
			}
		}
	}
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

bool LAMG::addAggregationLevel(Matrix &matrix, LAMGHierarchy &hierarchy, const count numTestVectors) {
	count iMax = 2;
	index i = 0;
	Vector B(iMax, std::numeric_limits<double>::max());
	std::vector<std::vector<index>> S(iMax, std::vector<index>(matrix.numberOfRows()));
	std::vector<count> numCoarseNodes(iMax);
	double alpha = 1.0;
	double delta = 0.9;
	double alpha_max = guard / cycleIndex;
	std::vector<Vector> testVectors = computeTestVectors(matrix, numTestVectors);
	Matrix C = computeAffinityMatrix(matrix, testVectors);

	std::vector<int64_t> status(matrix.numberOfRows(), -1);
	std::vector<count> aggregateSize(matrix.numberOfRows(), 1);

	addHighDegreeSeedNodes(matrix, status);

	while (alpha >= alpha_max && i < iMax) {
		delta *= 0.6;
		count nC;

		// aggregation stage
		aggregationStage(matrix, nC, C, testVectors, status, aggregateSize, delta);

		numCoarseNodes[i] = nC;
		alpha = nC / matrix.numberOfRows();
		alpha <= alpha_max ? B[i] = 1 - alpha :B[i] = 1 + alpha;

		std::unordered_map<index, index> old_to_new_ids;
		index newId = 0;
		for (index j = 0; j < matrix.numberOfRows(); ++j) {
			if (status[j] <= 0) {
				S[i][j] = newId;
				old_to_new_ids.insert({j, newId});
				newId++;
			} else {
				if (old_to_new_ids.find(j) != old_to_new_ids.end()) {
					S[i][j] = old_to_new_ids[j];
				} else {
					S[i][j] = newId;
					old_to_new_ids.insert({j, newId});
					newId++;
				}
			}
		}

		i++;
	}

	double min = std::numeric_limits<double>::max();
	for (index j = 0; j < iMax; ++j) {
		if (B[j] < min) {
			i = j;
			min = B[j];
		}
	}

	// create coarsened laplacian and interpolation matrix
	Matrix P(matrix.numberOfRows(), numCoarseNodes[i]);
	for (index j = 0; j < matrix.numberOfRows(); ++j) {
		P.setValue(j, S[i][j], 1);
	}

	Matrix coarseMatrix = Matrix::mTmMultiply(P, matrix * P);

	hierarchy.addLevel(matrix, P, LAMGHierarchy::AGGREGATION);

	return true;
}

LAMGHierarchy LAMG::buildHierarchy(const Graph &graph) {
	return LAMGHierarchy();
}

LAMGHierarchy LAMG::buildHierarchy(const Matrix &matrix) {
	return LAMGHierarchy();
}

} /* namespace NetworKit */
