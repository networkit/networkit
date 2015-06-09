///*
// * LAMG.cpp
// *
// *  Created on: 12.11.2014
// *      Author: Michael Wegner (michael.wegner@student.kit.edu)
// */
//
//#include "LAMG.h"
//
//namespace NetworKit {
//
//LAMG::LAMG(double guard, double cycleIndex) : guard(guard), cycleIndex(cycleIndex) {
//}
//
//count LAMG::lowDegreeNodes(const Matrix &matrix, std::vector<bool> &eliminate) const {
//	eliminate.resize(matrix.numberOfRows(), true);
//	count numEliminatedNodes = 0;
//	for (index i = 0; i < matrix.numberOfRows(); ++i) {
//		if (matrix.nnzInRow(i) - 1 <= 4 && eliminate[i]) { // node i has degree <= 4 and can be eliminated
//			numEliminatedNodes++;
//			matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight w){ // to maintain independence, mark all neighbors as not eliminated
//				if (j != i)	eliminate[j] = false;
//			});
//		} else {
//			eliminate[i] = false;
//		}
//	}
//
//	return numEliminatedNodes;
//}
//
//std::vector<Vector> LAMG::computeTestVectors(const Matrix &matrix, const count numOfVectors) const {
//	std::vector<Vector> testVectors(numOfVectors, Vector(matrix.numberOfColumns()));
//	Vector zeroVector(matrix.numberOfColumns(), 0.0);
//	for (count i = 0; i < numOfVectors; ++i) {
//		for (count j = 0; j < matrix.numberOfColumns(); ++j) {
//			testVectors[i][j] = Aux::Random::real(-1.0, 1.0);
//		}
//
//		// do 4 GS sweeps on the system Av = 0
//		testVectors[i] = smoother.relax(matrix, zeroVector, testVectors[i], 4);
//	}
//
//	return testVectors;
//}
//
//bool LAMG::fastRelaxationSpeed(const Matrix &matrix) const {
//	// check relaxation speed
//	Vector zeroVector(matrix.numberOfRows(), 0.0);
//	Vector x(matrix.numberOfRows());
////#pragma omp parallel for
//	for (index i = 0; i < matrix.numberOfRows(); ++i) {
//		x[i] = Aux::Random::real(-1,1);
//	}
//
//	// do 8 GS sweeps
//	x = smoother.relax(matrix, zeroVector, x, 8);
//
//	// do one more GS sweep
//	Vector y = smoother.relax(matrix, zeroVector, x, 1);
//
//	return (y - y.mean()).length() / (x - x.mean()).length() <= 0.7;
//}
//
//// TODO: optimize!
//Matrix LAMG::computeAffinityMatrix(const Matrix &matrix, const std::vector<Vector> &testVectors) const {
//	assert(testVectors.size() > 0);
//	Matrix C(matrix.numberOfRows(), matrix.numberOfColumns());
//	matrix.forNonZeroElementsInRowOrder([&](index i, index j, edgeweight value){
//		if (j >= i) {  // symmetric
//			double ij = 0.0;
//			double ii = 0.0;
//			double jj = 0.0;
//			for (count k = 0; k < testVectors.size(); ++k) {
//				ii += testVectors[k][i] * testVectors[k][i];
//				jj += testVectors[k][j] * testVectors[k][j];
//				ij += testVectors[k][i] * testVectors[k][j];
//			}
//
//			double weight = (ij * ij) / (ii * jj);
//			if (!(0.0 <= weight && weight <= 1.0)) {
//				std::runtime_error("weight out of range");
//			}
//
//
//			C.setValue(i, j, weight);
//			if (j > i) {
//				C.setValue(j, i, weight); // symmetric
//			}
//		}
//	});
//
//	return C;
//}
//
//void LAMG::addHighDegreeSeedNodes(const Matrix &matrix, std::vector<int64_t> &status) const {
//	std::vector<double> degree(matrix.numberOfRows(), 0.0);
////#pragma omp parallel for
//	for (index i = 0; i < matrix.numberOfRows(); ++i) {
//		double sum = 0.0;
//		matrix.forNonZeroElementsInRow(i, [&](index i, index j, double value){
//			if (j != i) {
//				degree[i] += -value * (matrix.nnzInRow(j) - 1);
//				sum += -value;
//			}
//		});
//
//		degree[i] += (matrix(i,i) - sum);
//		degree[i] /= matrix(i,i);
//	}
//
//	for (index i = 0; i < matrix.numberOfRows(); ++i) {
//		if (matrix.nnzInRow(i) - 1 >= 8 * degree[i]) {
//			status[i] = i; // mark node as seed node
//		}
//	}
//}
//
//std::vector<std::vector<index>> LAMG::computeStrongNeighbors(const Matrix &affinityMatrix, const double delta, const std::vector<int64_t> &status, std::vector<std::vector<index>> &bins) const {
//	std::vector<std::vector<index>> strongNeighbors(affinityMatrix.numberOfRows());
//	std::vector<double> maxNeighbor(affinityMatrix.numberOfRows(), std::numeric_limits<double>::min());
//	double overallMax = 0.0;
//	double overallMin = std::numeric_limits<double>::max();
//
//	affinityMatrix.forNonZeroElementsInRowOrder([&](index i, index j, double value) { // determine the highest affinity neighbor of each node
//		if (j != i && value > maxNeighbor[i]) {
//			maxNeighbor[i] = value;
//		}
//	});
//
//
//	affinityMatrix.forNonZeroElementsInRowOrder([&](index i, index j, double value) {
//		if (j != i && value >= delta * std::max(maxNeighbor[i], maxNeighbor[j])) {
//			strongNeighbors[i].push_back(j);
//		} else if (j == i) {
//			if (maxNeighbor[i] > overallMax) {
//				overallMax = maxNeighbor[i];
//			}
//			if (maxNeighbor[i] < overallMin) {
//				overallMin = maxNeighbor[i];
//			}
//		}
//	});
//
//	double h = fabs(overallMax - overallMin) < 1e-15? 1.0 : (double) bins.size() / (overallMax - overallMin);
//	for (index i = 0; i < affinityMatrix.numberOfRows(); ++i) {
//		if (status[i] == -1 && strongNeighbors[i].size() > 0) { // undecided nodes with strong neighbors
//			index binIndex = (index) floor(h * (maxNeighbor[i] - overallMin));
//			if (binIndex == bins.size()) { // last interval is closed on the right
//				binIndex--;
//			}
//
//			bins[binIndex].push_back(i);
//		}
//	}
//
//	return strongNeighbors;
//}
//
//bool LAMG::findBestSeed(const Matrix &affinityMatrix, const std::vector<index> &strongNeighborsOfU, const std::vector<int64_t> &status, const index u, index &s) const {
//	double maxAffinity = std::numeric_limits<double>::min();
//	for (index i = 0; i < strongNeighborsOfU.size(); ++i) {
//		index v = strongNeighborsOfU[i];
//		if (status[v] < 0 || status[v] == v) { // neighbor is seed or undecided
//			if (affinityMatrix(u, v) > maxAffinity) {
//				s = v;
//				maxAffinity = affinityMatrix(u, v);
//			}
//		}
//	}
//
//	return maxAffinity != std::numeric_limits<double>::min(); // we have found a strong neighbor which is seed or undecided
//}
//
//bool LAMG::findBestSeedEnergyCorrected(const Matrix &matrix, const Matrix &affinityMatrix, const std::vector<index> &strongNeighborsOfU, const std::vector<Vector> &testVectors, const std::vector<int64_t> &status, const index u, index &s) const {
//	Vector a(testVectors.size());
//	Vector b(testVectors.size());
//	Vector e(testVectors.size());
//	double Auu = matrix(u,u);
//	double maxAffinity = std::numeric_limits<double>::min();
//
//
//	for (index k = 0; k < testVectors.size(); ++k) {
//		matrix.forNonZeroElementsInRow(u, [&](index u, index v, edgeweight w){
//			if (u != v) {
//				a[k] += -w * testVectors[k][v] * testVectors[k][v];
//				b[k] += -w * testVectors[k][v];
//			}
//		});
//
//		a[k] *= 0.5;
//		double y = b[k] / Auu;
//		e[k] = (0.5*Auu*y - b[k])*y + a[k];
//	}
//
//	double minAggregateSize = std::numeric_limits<double>::max();
//	for (index t : strongNeighborsOfU) {
//		if (status[t] < 0 || status[t] == t) { // neighbor is seed or undecided
//			double maxMu = -1.0;
//			double ratioMax = 2.5;
//			for (index k = 0; k < testVectors.size(); ++k) {
//				double Ec = (0.5*Auu*testVectors[k][t] - b[k]) * testVectors[k][t] + a[k];
//
//				double mu = Ec / (e[k] + 1e15);
//				if (mu > maxMu) {
//					maxMu = mu;
//				}
//
//				if (maxMu > ratioMax) {
//					break;
//				}
//			}
//
//			if (maxMu <= ratioMax) {
//				double affinity = affinityMatrix(u,t);
//				if (affinity > maxAffinity) {
//					s = t;
//					maxAffinity = affinity;
//				}
//			}
//		}
//	}
//
//
//	return maxAffinity != std::numeric_limits<double>::min(); // we have found a strong neighbor which is seed or undecided
//}
//
//void LAMG::aggregationStage(const Matrix &matrix, count &numCoarseNodes, const Matrix &affinityMatrix, std::vector<Vector> &testVectors, std::vector<int64_t> &status, std::vector<count> &aggregateSize, double delta) const {
//	std::vector<std::vector<index>> bins(10);
//	std::vector<std::vector<index>> strongNeighbors = computeStrongNeighbors(affinityMatrix, delta, status, bins);
//
//	for (index k = bins.size() - 1; k > 0; k--) { // iterate over undecided nodes with strong neighbors in decreasing order of strongest neighbor
//		for (index i : bins[k]) {
//			if (status[i] == -1) { // node is still undecided
//				index s = 0;
//				//if (findBestSeed(affinityMatrix, strongNeighbors[i], status, i, s)) {
//				if (findBestSeedEnergyCorrected(matrix, affinityMatrix, strongNeighbors[i], testVectors, status, i, s)) {
//					status[s] = s; // s becomes seed
//					status[i] = s; // i's seed is s
//					numCoarseNodes--;
//
//					for (index j = 0; j < testVectors.size(); ++j) { // update test vectors
//						testVectors[j][i] = testVectors[j][s];
//					}
//
//					aggregateSize[s]++;
//					aggregateSize[i] = aggregateSize[s];
//				}
//			}
//		}
//
//		if (numCoarseNodes <= matrix.numberOfRows() * guard / cycleIndex) {
//			break;
//		}
//	} // iterate over bins
//}
//
//bool LAMG::addEliminationLevel(Matrix &matrix, LAMGHierarchy &hierarchy) const {
//	INFO("Add Elimination level");
//	bool nodesEliminated = false;
//	count iterations = 0;
//
//	Matrix P = DiagonalMatrix(matrix.numberOfRows()); // initial interpolation matrix is unit matrix
//	Matrix Q(matrix.numberOfRows());
//	while (matrix.numberOfRows() > 1 && iterations < 1000) {
//		std::vector<bool> eliminate;
//		count numFNodes = lowDegreeNodes(matrix, eliminate); // find independent nodes with deg <= 4
//		if (numFNodes < 0.01 * matrix.numberOfRows() || matrix.numberOfRows() < 200) {
//			if (nodesEliminated) { // add level if nodes have been eliminated
//				hierarchy.addLevel(matrix, P, Q);
//			}
//			return nodesEliminated;
//		}
//
//		count numCNodes = matrix.numberOfRows() - numFNodes; // number of coarse nodes
//		assert(numCNodes > 0);
//
//		// calculate matrix S = {(-A_FC^T * A_FF^-1)^T, Ic}^T and permutation matrix M
//		Matrix S(matrix.numberOfRows(), numCNodes);
//		Matrix M(matrix.numberOfRows());
//		Matrix Qi(matrix.numberOfRows());
//
//		// order nodes s.t. F nodes appear first and then all other nodes (C nodes)
//		std::vector<index> sIndex(matrix.numberOfRows());
//		index fIndex = 0;
//		index cIndex = numFNodes;
//		for (index i = 0; i < matrix.numberOfRows(); ++i) {
//			if (eliminate[i]) { // node i is an F node
//				sIndex[i] = fIndex;
//				M.setValue(i, fIndex, 1);
//				Qi.setValue(i, fIndex, 1 / matrix(i,i));
//				fIndex++;
//			} else { // node i is a C node
//				sIndex[i] = cIndex;
//				M.setValue(i, cIndex, 1);
//				cIndex++;
//			}
//		}
//
//		//INFO("fIndex=", fIndex, ", cIndex=", cIndex);
//
//		for (index i = 0; i < matrix.numberOfRows(); ++i) {
//			if (eliminate[i]) { // node i is an F node
//				double weightedDegree = matrix(i,i);
//				matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight value){ // (-A_FC^T * A_FF^-1)^T in row i
//					if (j != i)	S.setValue(sIndex[i], sIndex[j] - numFNodes, -value / weightedDegree);
//					if (j != i) assert(!eliminate[j]);
//				});
//			} else { // node i is a C node
//				S.setValue(sIndex[i], sIndex[i] - numFNodes, 1);
//			}
//		}
//
//		// compute current interpolation matrix P_iteration
//		Matrix P_iteration = M * S;
//
//		// compute coarsened matrix
//		matrix = Matrix::mTmMultiply(P_iteration, matrix * P_iteration);
//
//		// check that sum of each row is zero (matrix is laplacian)
//		for (index i = 0; i < matrix.numberOfRows(); ++i) {
//			if (abs(matrix(i,i)) <= 1e-15) INFO("Zero entry in coarse elimination matrix!");
//			edgeweight sum = 0.0;
//			matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight value) {
//				sum += value;
//			});
//
//			matrix.setValue(i, i, matrix(i,i) - sum);
//		}
//
//		Q += Matrix::mmTMultiply(P * Qi, P);
//		P = P * P_iteration;
//		nodesEliminated = true;
//		iterations++;
//	}
//
//	if (nodesEliminated) { // add level if nodes have been eliminated
//		hierarchy.addLevel(matrix, P, Q);
//	}
//
//	return nodesEliminated;
//}
//
//bool LAMG::addAggregationLevel(Matrix &matrix, LAMGHierarchy &hierarchy, const count numTestVectors) const {
//	count iMax = 2;
//	index i = 0;
//	Vector B(iMax, std::numeric_limits<double>::max());
//	std::vector<std::vector<index>> S(iMax, std::vector<index>(matrix.numberOfRows(), std::numeric_limits<index>::max()));
//	std::vector<count> numCoarseNodes(iMax, matrix.numberOfRows());
//	double alpha = 1.0;
//	double delta = 0.9;
//	double alpha_max = guard / cycleIndex;
//	std::vector<Vector> testVectors = computeTestVectors(matrix, numTestVectors);
//	Matrix affinityMatrix = computeAffinityMatrix(matrix, testVectors);
//
//	std::vector<int64_t> status(matrix.numberOfRows(), -1);
//	std::vector<count> aggregateSize(matrix.numberOfRows(), 1);
//
//	addHighDegreeSeedNodes(matrix, status);
//
//	while (i < 1 || (alpha >= alpha_max && i < iMax)) {
//		count nC = i > 0? numCoarseNodes[i - 1] : matrix.numberOfRows();
//
//		// aggregation stage
//		aggregationStage(matrix, nC, affinityMatrix, testVectors, status, aggregateSize, delta);
//
//		alpha = nC / (double) matrix.numberOfRows();
//		alpha <= alpha_max ? B[i] = 1 - alpha : B[i] = 1 + alpha;
//
//		index newId = 0;
//		count undecided = 0;
//		count seedNodes = 0;
//
//		count decided = 0;
//
//		for (index j = 0; j < matrix.numberOfRows(); ++j) {
//			if (status[j] < 0) { // undecided nodes which will remain in the coarse set
//				S[i][j] = newId;
//				newId++;
//				undecided++;
//			} else {
//				if (status[j] == j && S[i][j] > newId) { // seed node which has not been added yet
//					S[i][j] = newId;
//					newId++;
//					seedNodes++;
//				} else if (status[j] != j) { // aggregated node
//					if (S[i][status[j]] > newId) { // seed of j has not been added to S yet, therefore add it now
//						S[i][status[j]] = newId;
//						newId++;
//						seedNodes++;
//					}
//
//					S[i][j] = S[i][status[j]];
//					decided++;
//				}
//			}
//		}
//
//		//INFO("added aggregate set S", i);
//		numCoarseNodes[i] = nC;
//		i++;
//		delta *= 0.6;
//	}
//
//	double min = std::numeric_limits<double>::max();
//	for (index j = 0; j < iMax; ++j) {
//		if (B[j] < min) {
//			i = j;
//			min = B[j];
//		}
//	}
//
//	if (numCoarseNodes[i] < matrix.numberOfRows()) {
//		// create interpolation matrix
//		Matrix P(matrix.numberOfRows(), numCoarseNodes[i]);
//		for (index j = 0; j < matrix.numberOfRows(); ++j) {
//			P.setValue(j, S[i][j], 1);
//		}
//
//		// create coarsened laplacian
//		matrix = Matrix::mTmMultiply(P, matrix * P);
//
//		// check that sum of each row is zero (matrix is laplacian)
//		for (index i = 0; i < matrix.numberOfRows(); ++i) {
//			edgeweight sum = 0.0;
//			matrix.forNonZeroElementsInRow(i, [&](index i, index j, edgeweight value){
//				sum += value;
//			});
//
//			sum -= matrix(i,i);
//			matrix.setValue(i, i, -sum);
//		}
//		hierarchy.addLevel(matrix, P, LAMGHierarchy::AGGREGATION);
//		return true;
//	}
//
//	return false;
//}
//
//LAMGHierarchy LAMG::buildHierarchy(const Matrix &matrix) const {
//	Matrix A = matrix;
//	LAMGHierarchy hierarchy(A);
//	count numAggregations = 0;
//	INFO("FINEST\t", matrix.numberOfRows(), "\t", matrix.nnz() - matrix.numberOfRows());
//	while (!fastRelaxationSpeed(A) && A.numberOfRows() > 200) {
//		// first search for low degree nodes to eliminate
//		if (addEliminationLevel(A, hierarchy)) {
//			INFO(hierarchy.getNumLevels() - 1, " ELIM\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
//			if (fastRelaxationSpeed(A) || A.numberOfRows() <= 200) break;
//		}
//
//		// aggregation
//		if (addAggregationLevel(A, hierarchy, 4 + numAggregations)) {
//			INFO(hierarchy.getNumLevels() - 1, " AGG\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
//			numAggregations++;
//		}
//	}
//
//	return hierarchy;
//}
//
//} /* namespace NetworKit */
