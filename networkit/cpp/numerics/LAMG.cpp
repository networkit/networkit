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

count LAMG::lowDegreeNodes(const Matrix &matrix, std::vector<bool> &eliminate) const {
	eliminate.resize(matrix.numberOfRows(), true);
	count numEliminatedNodes = 0;
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (matrix.nnzInRow(i) - 1 <= 4 && eliminate[i]) { // node i has degree <= 4 and can be eliminated
			numEliminatedNodes++;
			matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight w){ // to maintain independence, mark all neighbors as not eliminated
				if (j != i)	eliminate[j] = false;
			});
		} else {
			eliminate[i] = false;
		}
	}

	return numEliminatedNodes;
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

bool LAMG::fastRelaxationSpeed(const Matrix &matrix) const {
	// check relaxation speed
	Vector zeroVector(matrix.numberOfRows(), 0.0);
	Vector testVector(matrix.numberOfRows());
//#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		testVector[i] = Aux::Random::real(-1,1);
	}

	// do 14 GS sweeps
	testVector = smoother.relax(matrix, zeroVector, testVector, 14);

	// do one more GS sweep
	Vector iterTestVector = smoother.relax(matrix, zeroVector, testVector, 1);

	return iterTestVector.length() / testVector.length() <= 0.7;
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

			double weight = (ij * ij) / (ii * jj);

			C.setValue(i, j, weight);
			if (j > i) {
				C.setValue(j, i, weight); // symmetric
			}
		}
	});

	return C;
}

void LAMG::addHighDegreeSeedNodes(const Matrix &matrix, std::vector<int64_t> &status) const {
	std::vector<count> degree(matrix.numberOfRows(), 0);
//#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		degree[i] = matrix.nnzInRow(i) - 1;
	}

	std::sort(degree.begin(), degree.end());
	count median = degree[floor(matrix.numberOfRows() / 2)];

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (matrix.nnzInRow(i) - 1 >= 8 * median) {
			status[i] = i;
		}
	}
}

std::vector<std::vector<index>> LAMG::computeStrongNeighbors(const Matrix &affinityMatrix, const double delta) const {
	std::vector<std::vector<index>> strongNeighbors(affinityMatrix.numberOfRows());
	std::vector<double> maxNeighbor(affinityMatrix.numberOfRows(), std::numeric_limits<double>::min());

	affinityMatrix.forNonZeroElementsInRowOrder([&](index i, index j, double value) { // determine the highest affinity neighbor of each node
		if (j != i && value > maxNeighbor[i]) {
			maxNeighbor[i] = value;
		}
	});


	affinityMatrix.forNonZeroElementsInRowOrder([&](index i, index j, double value) {
		if (j != i && value >= delta * std::max(maxNeighbor[i], maxNeighbor[j])) {
			strongNeighbors[i].push_back(j);
		}
	});

	return strongNeighbors;
}

bool LAMG::findBestSeed(const Matrix &affinityMatrix, const std::vector<index> &strongNeighborsOfU, const std::vector<int64_t> &status, const index u, index &s) const {
	double maxAffinity = std::numeric_limits<double>::min();
	for (index i = 0; i < strongNeighborsOfU.size(); ++i) {
		index v = strongNeighborsOfU[i];
		if (status[v] < 0 || status[v] == v) { // neighbor is seed or undecided
			if (affinityMatrix(u, v) > maxAffinity) {
				s = v;
				maxAffinity = affinityMatrix(u, v);
			}
		}
	}

	return maxAffinity != std::numeric_limits<double>::min(); // we have found a strong neighbor which is seed or undecided
}

void LAMG::aggregationStage(const Matrix &matrix, count &numCoarseNodes, const Matrix &affinityMatrix, std::vector<Vector> &testVectors, std::vector<int64_t> &status, std::vector<count> &aggregateSize, double delta) const {
	std::vector<std::vector<index>> strongNeighbors = computeStrongNeighbors(affinityMatrix, delta);
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (status[i] == -1 && strongNeighbors[i].size() > 0) { // undecided nodes with strong neighbors
			index s;
			if (findBestSeed(affinityMatrix, strongNeighbors[i], status, i, s)) {
				status[s] = s; // s becomes seed
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

bool LAMG::addEliminationLevel(Matrix &matrix, LAMGHierarchy &hierarchy) const {
	bool nodesEliminated = false;

	Matrix P = DiagonalMatrix(matrix.numberOfRows()); // initial interpolation matrix is unit matrix
	while (matrix.numberOfRows() > 1) {
		std::vector<bool> eliminate;
		count numFNodes = lowDegreeNodes(matrix, eliminate); // find independent nodes with deg <= 4
		if (numFNodes < 0.01 * matrix.numberOfRows()) {
			if (nodesEliminated) { // add level if nodes have been eliminated
				hierarchy.addLevel(matrix, P, LAMGHierarchy::ELIMINATION);
			}
			return nodesEliminated;
		}

		count numCNodes = matrix.numberOfRows() - numFNodes; // number of coarse nodes

		INFO("Elimination-Stats: numFNodes=", numFNodes, " numCNodes=", numCNodes);

		// calculate matrix S = {(-A_FC^T * A_FF^-1)^T, Ic} and permutation matrix M
		Matrix S(matrix.numberOfRows(), numCNodes);
		Matrix M(matrix.numberOfRows());

		// order nodes s.t. F nodes appear first and then all other nodes (C nodes)
		std::vector<index> sIndex(matrix.numberOfRows());
		index fIndex = 0;
		index cIndex = numFNodes;
		for (index i = 0; i < matrix.numberOfRows(); ++i) {
			if (eliminate[i]) {
				sIndex[i] = fIndex;
				M.setValue(i, fIndex, 1);
				fIndex++;
			} else {
				sIndex[i] = cIndex;
				M.setValue(i, cIndex, 1);
				cIndex++;
			}
		}

		INFO("fIndex=", fIndex, ", cIndex=", cIndex);

		for (index i = 0; i < matrix.numberOfRows(); ++i) {
			if (eliminate[i]) {
				matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight value){ // (-A_FC^T * A_FF^-1)^T in row i
					if (j != i)	S.setValue(sIndex[i], sIndex[j] - numFNodes, -value / matrix(i,i));
				});
			} else {
				S.setValue(sIndex[i], sIndex[i] - numFNodes, 1);
			}
		}

		// compute current interpolation matrix P_iteration
		Matrix P_iteration = M * S;

		// compute coarsened matrix
		matrix = Matrix::mTmMultiply(P_iteration, matrix * P_iteration);
		P = P * P_iteration;
		nodesEliminated = true;
	}

	if (nodesEliminated) { // add level if nodes have been eliminated
		hierarchy.addLevel(matrix, P, LAMGHierarchy::ELIMINATION);
	}

	return nodesEliminated;
}

bool LAMG::addAggregationLevel(Matrix &matrix, LAMGHierarchy &hierarchy, const count numTestVectors) const {
	count iMax = 2;
	index i = 0;
	Vector B(iMax, std::numeric_limits<double>::max());
	std::vector<std::vector<index>> S(iMax, std::vector<index>(matrix.numberOfRows(), std::numeric_limits<index>::max()));
	std::vector<count> numCoarseNodes(iMax, matrix.numberOfRows());
	double alpha = 1.0;
	double delta = 0.9;
	double alpha_max = guard / cycleIndex;
	std::vector<Vector> testVectors = computeTestVectors(matrix, numTestVectors);
	Matrix C = computeAffinityMatrix(matrix, testVectors);

	std::vector<int64_t> status(matrix.numberOfRows(), -1);
	std::vector<count> aggregateSize(matrix.numberOfRows(), 1);

	addHighDegreeSeedNodes(matrix, status);

	while (alpha >= alpha_max && i < iMax) {
		INFO("Calculating aggregate set S", i);
		count nC = i > 0? numCoarseNodes[i - 1] : matrix.numberOfRows();
		INFO("initial coarse nodes: ", nC);

		// aggregation stage
		aggregationStage(matrix, nC, C, testVectors, status, aggregateSize, delta);

		alpha = nC / (double) matrix.numberOfRows();
		alpha <= alpha_max ? B[i] = 1 - alpha : B[i] = 1 + alpha;
		INFO("alpha=", alpha);

		index newId = 0;
		count undecided = 0;
		count seedNodes = 0;

		count decided = 0;

		for (index j = 0; j < matrix.numberOfRows(); ++j) {
			if (status[j] < 0) { // undecided nodes
				S[i][j] = newId;
				newId++;
				undecided++;
			} else {
				if (status[j] == j && S[i][j] > newId) { // seed node which has not been added yet
					S[i][j] = newId;
					newId++;
					seedNodes++;
				} else if (status[j] != j) { // aggregated node
					if (S[i][status[j]] > newId) { // seed of j has not been added to S yet, therefore add it now
						S[i][status[j]] = newId;
						newId++;
						seedNodes++;
					}

					S[i][j] = S[i][status[j]];
					decided++;
				}
			}
		}

		INFO("numCoarseNodes=", nC, " and newId=", newId, " undecided=", undecided, " seeds=", seedNodes, " decided=", decided);

		INFO("added aggregate set S", i);
		numCoarseNodes[i] = nC;
		i++;
		delta *= 0.6;
	}

	INFO("Look for best aggregate");
	double min = std::numeric_limits<double>::max();
	for (index j = 0; j < iMax; ++j) {
		if (B[j] < min) {
			i = j;
			min = B[j];
		}
	}

	if (numCoarseNodes[i] < matrix.numberOfRows()) {
		INFO("Create interpolation matrix");
		// create coarsened laplacian and interpolation matrix
		Matrix P(matrix.numberOfRows(), numCoarseNodes[i]);
		for (index j = 0; j < matrix.numberOfRows(); ++j) {
			P.setValue(j, S[i][j], 1);
		}

		INFO("Create coarse laplacian");
		matrix = Matrix::mTmMultiply(P, matrix * P);

		INFO("Adding hierarchy level");
		hierarchy.addLevel(matrix, P, LAMGHierarchy::AGGREGATION);

		INFO("Finished aggregation level");

		return true;
	}

	return false;
}

MultigridHierarchy LAMG::buildHierarchy(const Matrix &matrix) const {
	Matrix A = matrix;
	LAMGHierarchy hierarchy(A);
	count numAggregations = 0;
	INFO("LAMG-Setup");
	while (!fastRelaxationSpeed(A) && A.numberOfRows() > 150) {
		// first search for low degree nodes to eliminate
		INFO("Searching for nodes to eliminate on level ", hierarchy.getNumLevels());
		if (addEliminationLevel(A, hierarchy)) {
			INFO("Elimination level added");
		}

		// aggregation
		INFO("Aggregating nodes on level ", hierarchy.getNumLevels());
		if (addAggregationLevel(A, hierarchy, 8 + numAggregations)) {
			numAggregations++;
			INFO("Aggregation level added");
		}
	}
	INFO("DONE");

	return hierarchy;
}

} /* namespace NetworKit */
