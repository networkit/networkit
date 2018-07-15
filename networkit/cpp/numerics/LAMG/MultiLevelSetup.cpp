/*
 * MultiLevelSetup.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultiLevelSetup.h"

namespace NetworKit {

template<>
void MultiLevelSetup<CSRMatrix>::setup(const CSRMatrix& matrix, LevelHierarchy<CSRMatrix>& hierarchy) const {
	CSRMatrix A = matrix;
	A.sort();
	setupForMatrix(A, hierarchy);
}

template<>
void MultiLevelSetup<CSRMatrix>::computeStrongAdjacencyMatrix(const CSRMatrix& matrix, CSRMatrix& strongAdjMatrix) const {
	std::vector<double> maxNeighbor(matrix.numberOfRows(), std::numeric_limits<double>::min());
#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
		matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (i != j && -value > maxNeighbor[i]) {
				maxNeighbor[i] = -value;
			}
		});
	}

	std::vector<index> rowIdx(matrix.numberOfRows()+1, 0);
	matrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) {
		if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
			++rowIdx[i+1];
		}
	});

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		rowIdx[i+1] += rowIdx[i];
	}

	count nnz = rowIdx[matrix.numberOfRows()];
	std::vector<index> columnIdx(nnz);
	std::vector<double> nonZeros(nnz);

#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
		index cIdx = rowIdx[i];
		matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
				columnIdx[cIdx] = j;
				nonZeros[cIdx] = -value;
				++cIdx;
			}
		});
	}

	strongAdjMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), rowIdx, columnIdx, nonZeros, 0.0, matrix.sorted());
}

template<>
void MultiLevelSetup<CSRMatrix>::computeAffinityMatrix(const CSRMatrix& matrix, const std::vector<Vector> &tVs, CSRMatrix& affinityMatrix) const {
	assert(tVs.size() > 0);

	std::vector<index> rowIdx(matrix.numberOfRows()+1);
	std::vector<index> columnIdx(matrix.nnz());
	std::vector<double> nonZeros(matrix.nnz());

#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
		rowIdx[i+1] = matrix.nnzInRow(i);
	}

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		rowIdx[i+1] += rowIdx[i];
	}

	std::vector<double> normSquared(matrix.numberOfRows(), 0.0);
#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
		for (index k = 0; k < tVs.size(); ++k) {
			normSquared[i] += tVs[k][i] * tVs[k][i];
		}
	}

#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(matrix.numberOfRows()); ++i) {
		double nir = 1.0 / normSquared[i];
		index cIdx = rowIdx[i];
		matrix.forNonZeroElementsInRow(i, [&](index j, double val) {
			double ij = 0.0;
			for (index k = 0; k < tVs.size(); ++k) {
				ij += tVs[k][i] * tVs[k][j];
			}

			double value = (ij * ij) * nir / normSquared[j];
			columnIdx[cIdx] = j;
			nonZeros[cIdx] = value;
			++cIdx;
		});
	}

	affinityMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), rowIdx, columnIdx, nonZeros, 0.0, matrix.sorted());
}

template<>
void MultiLevelSetup<CSRMatrix>::galerkinOperator(const CSRMatrix& P, const CSRMatrix& A, const std::vector<index>& PColIndex, const std::vector<std::vector<index>>& PRowIndex, CSRMatrix& B) const {
	std::vector<Triplet> triplets;
	SparseAccumulator spa(P.numberOfColumns());
	for (index i = 0; i < P.numberOfColumns(); ++i) {
		for (index k : PRowIndex[i]) {
			double Pki = P(k,i);
			A.forNonZeroElementsInRow(k, [&](index l, double value) {
				index j = PColIndex[l];
				spa.scatter(Pki * value * P(l, j), j);
			});
		}

		spa.gather([&](index i, index j, double value) {
			triplets.push_back({i,j,value});
		});

		spa.increaseRow();
	}

	B = CSRMatrix(P.numberOfColumns(), P.numberOfColumns(), triplets, 0.0, true);
}


template<>
void MultiLevelSetup<CSRMatrix>::eliminationOperators(const CSRMatrix& matrix, const std::vector<index>& fSet, const std::vector<index>& coarseIndex, CSRMatrix& P, Vector& q) const {
	std::vector<Triplet> triples;
	q = Vector(fSet.size());
	for (index k = 0; k < fSet.size(); ++k) { // Afc
		matrix.forNonZeroElementsInRow(fSet[k], [&](index j, edgeweight w){
			if (fSet[k] == j) {
				q[k] = 1.0 / w;
			} else {
				triples.push_back({k, coarseIndex[j], w});
			}
		});
	}

	for (index i = 0; i < triples.size(); ++i) { // * -Aff^-1
		triples[i].value *= -q[triples[i].row];
	}

	P = CSRMatrix(fSet.size(), coarseIndex.size() - fSet.size(), triples, 0.0, matrix.sorted());
}

template<>
void MultiLevelSetup<CSRMatrix>::coarseningAggregation(CSRMatrix& matrix, LevelHierarchy<CSRMatrix>& hierarchy, Vector& tv, count numTVVectors) const {
	Vector B(SETUP_MAX_AGGREGATION_STAGES, std::numeric_limits<double>::max());
	std::vector<std::vector<index>> S(SETUP_MAX_AGGREGATION_STAGES, std::vector<index>(matrix.numberOfRows(), std::numeric_limits<index>::max()));
	std::vector<index> status(matrix.numberOfRows(), UNDECIDED);
	std::vector<count> nc(SETUP_MAX_AGGREGATION_STAGES, matrix.numberOfRows());

	double alpha = 1.0;
	double maxCoarseningRatio = SETUP_COARSENING_WORK_GUARD / SETUP_CYCLE_INDEX;
	count stage = 0;
	count nC = matrix.numberOfRows();


	// generate TVs
	std::vector<Vector> tVs = generateTVs(matrix, tv, numTVVectors);

	// compute strong adjacency matrix
	CSRMatrix Wstrong;
	computeStrongAdjacencyMatrix(matrix, Wstrong);

	// compute affinityMatrix
	CSRMatrix affinityMatrix;
	computeAffinityMatrix(Wstrong, tVs, affinityMatrix);

	// mark all locally high-degree nodes as seeds
	addHighDegreeSeedNodes(matrix, status);

	// aggregate all loose nodes
	aggregateLooseNodes(Wstrong, status, nC);

	nc[0] = nC;
	while (stage < SETUP_MIN_AGGREGATION_STAGES || (alpha >= maxCoarseningRatio && stage < SETUP_MAX_AGGREGATION_STAGES)) {
		nC = stage > 0? nc[stage - 1] : nc[0];

		// aggregation stage
		aggregationStage(matrix, nC, Wstrong, affinityMatrix, tVs, status);

		alpha = (double) nC / (double) matrix.numberOfRows();
		alpha <= maxCoarseningRatio? B[stage] = 1.0-alpha : B[stage] = 1.0+alpha;

		S[stage] = status;
		nc[stage] = nC;
		stage++;
	}

	double min = B[0];
	index bestAggregate = 0;
	for (index i = 1; i < stage; ++i) {
		if (B[i] < min) {
			bestAggregate = i;
			min = B[i];
		}
	}

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (S[bestAggregate][i] == UNDECIDED) { // undediced nodes become their own seeds
			S[bestAggregate][i] = i;
		}
	}

	std::vector<index> indexFine(matrix.numberOfRows(), 0);
	index newIndex = 0;
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if (S[bestAggregate][i] == i) {
			indexFine[i] = newIndex++;
		}
	}

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		status[i] = indexFine[S[bestAggregate][i]];
	}

	assert(newIndex == nc[bestAggregate]);

	// create interpolation matrix
	std::vector<Triplet> pTriples(matrix.numberOfRows());
	std::vector<Triplet> rTriples(matrix.numberOfRows());
	std::vector<index> PColIndex(matrix.numberOfRows());
	std::vector<std::vector<index>> PRowIndex(nc[bestAggregate]);

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		pTriples[i] = {i, status[i], 1};
		rTriples[i] = {status[i], i, 1};
		PColIndex[i] = status[i];
		PRowIndex[status[i]].push_back(i);
	}

	CSRMatrix P(matrix.numberOfRows(), nc[bestAggregate], pTriples, 0.0, matrix.sorted());
	CSRMatrix R(nc[bestAggregate], matrix.numberOfRows(), rTriples, 0.0, matrix.sorted());

	// create coarsened laplacian
	galerkinOperator(P, matrix, PColIndex, PRowIndex, matrix);

	hierarchy.addAggregationLevel(matrix, P, R);
}

} /* namespace NetworKit */
