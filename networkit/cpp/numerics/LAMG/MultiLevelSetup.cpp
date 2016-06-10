/*
 * MultiLevelSetup.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael Wegner (michael,wegner@student.kit.edu)
 */

#include "MultiLevelSetup.h"
#include "Level/LevelElimination.h"
#include "Level/EliminationStage.h"
#include "LAMGSettings.h"
#include "../../auxiliary/StringTools.h"
#include "../../auxiliary/Enforce.h"
#include "../../auxiliary/Timer.h"
#include "../../algebraic/CSRMatrix.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <set>
#include "omp.h"

namespace NetworKit {

#ifndef NDEBUG
count MultiLevelSetup::eliminationTime = 0;
count MultiLevelSetup::schurComplementTime = 0;
count MultiLevelSetup::aggregationTime = 0;
#endif

MultiLevelSetup::MultiLevelSetup(const Smoother &smoother) : smoother(smoother) {
}

void MultiLevelSetup::setup(const Graph &G, LevelHierarchy &hierarchy) const {
	setup(CSRMatrix::graphLaplacian(G), hierarchy);
}

void MultiLevelSetup::setup(const CSRMatrix &matrix, LevelHierarchy &hierarchy) const {
	CSRMatrix A = matrix;
	hierarchy.addFinestLevel(A);
#ifndef NDEBUG
	DEBUG("FINEST\t", matrix.numberOfRows(), "\t", matrix.nnz());
#endif

	bool doneCoarsening = false;
	count numTVs = TV_NUM;
	index level = 0;
	A.sort();
	while (!doneCoarsening) {
		// ELIMINATION
		if (coarseningElimination(A, hierarchy)) {
			if (!canCoarsen(A)) doneCoarsening = true;
			level++;
#ifndef NDEBUG
			DEBUG(level, " ELIM\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
#endif
		}

		// AGGREGATION
		Vector tv;
		if (doneCoarsening || isRelaxationFast(A, level, tv)) {
			doneCoarsening = true;
		} else {
			coarseningAggregation(A, hierarchy, tv, numTVs);
			level++;
#ifndef NDEBUG
			DEBUG(level, " AGG\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
#endif
			if (numTVs < TV_MAX) {
				numTVs += TV_INC;
			}
		}

		if (!canCoarsen(A)) doneCoarsening = true;
	}

	hierarchy.setLastAsCoarsest();

#ifndef NDEBUG
	DEBUG("Elimination: ", eliminationTime);
	DEBUG("Schur: ", schurComplementTime);
	DEBUG("Aggregation: ", aggregationTime);
#endif
}

bool MultiLevelSetup::coarseningElimination(CSRMatrix &matrix, LevelHierarchy &hierarchy) const {
#ifndef NDEBUG
	Aux::Timer elimTimer;
	Aux::Timer schurTimer;
	elimTimer.start();
#endif
	std::vector<EliminationStage> coarseningStages;
	count stageNum = 0;
	while (stageNum < SETUP_ELIMINATION_MAX_STAGES) {
		if (matrix.numberOfRows() <= MAX_DIRECT_SOLVE_SIZE) break; // we do not need to coarsen the matrix any further
		std::vector<bool> fNode;
		count nf = lowDegreeSweep(matrix, fNode, stageNum);
		count nc = matrix.numberOfRows() - nf;

		if (nc == 0) { // do not eliminate all nodes -> leave one entry in c
			nc = 1;
			nf--;
		}

		// add f nodes to fSet and c nodes to cSet
		std::vector<index> fSet(nf);
		std::vector<index> cSet(nc);

		std::vector<index> coarseIndex(matrix.numberOfRows());
		count numFNodes = 0;
		for (index i = 0, fIndex = 0, cIndex = 0; i < matrix.numberOfRows(); ++i) {
			if (fNode[i] && fIndex < nf) {
				coarseIndex[i] = fIndex;
				fSet[fIndex++] = i;
				numFNodes++;
			} else {
				coarseIndex[i] = cIndex;
				cSet[cIndex++] = i;
			}
		}

		if (nf <= SETUP_ELIMINATION_MIN_ELIM_FRACTION * matrix.numberOfRows()) {
			break;
		}

		CSRMatrix P;
		Vector q;
		eliminationOperators(matrix, fSet, coarseIndex, P, q);
		coarseningStages.push_back(EliminationStage(P, q, fSet, cSet));

#ifndef NDEBUG
		schurTimer.start();
#endif

		CSRMatrix Acc = matrix.subMatrix(cSet, cSet); // Schur complement
		CSRMatrix Acf = matrix.subMatrix(cSet, fSet); // Schur complement

		matrix = Acc + Acf * P;

#ifndef NDEBUG
		schurTimer.stop();
		schurComplementTime += schurTimer.elapsedMilliseconds();
#endif

		stageNum++;
	}

	if (stageNum != 0) { // we have coarsened the matrix
		hierarchy.addEliminationLevel(matrix, coarseningStages);
#ifndef NDEBUG
		elimTimer.stop();
		eliminationTime += elimTimer.elapsedMilliseconds();
		//schurComplementTime += schurTimer.elapsedMilliseconds();
#endif
		return true;
	}
#ifndef NDEBUG
	elimTimer.stop();
	eliminationTime += elimTimer.elapsedMilliseconds();
	//schurComplementTime += schurTimer.elapsedMilliseconds();
#endif

	return false;
}

count MultiLevelSetup::lowDegreeSweep(const CSRMatrix &matrix, std::vector<bool> &fNode, index stage) const {
	fNode.resize(matrix.numberOfRows(), true); // first mark all nodes as f nodes
	count numFNodes = 0;
	int degreeOffset = stage != 0;

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if ((int) matrix.nnzInRow(i) - degreeOffset <= (int)SETUP_ELIMINATION_MAX_DEGREE && fNode[i]) { // node i has degree <= 4 and can be eliminated
			numFNodes++;
			matrix.forNonZeroElementsInRow(i, [&](index j, edgeweight w){ // to maintain independence, mark all neighbors as not eliminated
				if (j != i)	{ // all neighbors of this f node are c nodes
					fNode[j] = false;
				}
			});
		} else { // node has high degree, thus it is a c node
			fNode[i] = false;
		}
	}

	return numFNodes;
}

void MultiLevelSetup::eliminationOperators(const CSRMatrix &matrix, const std::vector<index> &fSet, const std::vector<index> &coarseIndex, CSRMatrix &P, Vector &q) const {
	std::vector<CSRMatrix::Triple> triples;
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

	P = CSRMatrix(fSet.size(), coarseIndex.size() - fSet.size(), triples, matrix.sorted());
}

void MultiLevelSetup::subMatrix(const CSRMatrix &matrix, const std::vector<index> &rows, const std::vector<index> &columns, const std::vector<index> &coarseIndex, CSRMatrix &result) const {
	std::vector<CSRMatrix::Triple> triples;

	for (index k = 0; k < rows.size(); ++k) {
		matrix.forNonZeroElementsInRow(rows[k], [&](index j, edgeweight value) {
			if (coarseIndex[j] < columns.size() && columns[coarseIndex[j]] == j) { // check if neighbor is in columns
				triples.push_back({k, coarseIndex[j], value});
			}
		});
	}

	result = CSRMatrix(rows.size(), columns.size(), triples, matrix.sorted());
}

void MultiLevelSetup::coarseningAggregation(CSRMatrix &matrix, LevelHierarchy &hierarchy, Vector &tv, count numTVVectors) const {
#ifndef NDEBUG
	Aux::Timer aggTimer;
	aggTimer.start();
#endif
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
	std::vector<CSRMatrix::Triple> pTriples(matrix.numberOfRows());
	std::vector<CSRMatrix::Triple> rTriples(matrix.numberOfRows());
	std::vector<index> PColIndex(matrix.numberOfRows());
	std::vector<std::vector<index>> PRowIndex(nc[bestAggregate]);

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		pTriples[i] = {i, status[i], 1};
		rTriples[i] = {status[i], i, 1};
		PColIndex[i] = status[i];
		PRowIndex[status[i]].push_back(i);
	}

	CSRMatrix P(matrix.numberOfRows(), nc[bestAggregate], pTriples, matrix.sorted());
	CSRMatrix R(nc[bestAggregate], matrix.numberOfRows(), rTriples, matrix.sorted());

	// create coarsened laplacian
	galerkinOperator(P, matrix, PColIndex, PRowIndex, matrix);

	hierarchy.addAggregationLevel(matrix, P, R);

#ifndef NDEBUG
	aggTimer.stop();
	aggregationTime += aggTimer.elapsedMilliseconds();
#endif
}

std::vector<Vector> MultiLevelSetup::generateTVs(const CSRMatrix &matrix, Vector &tv, count numVectors) const {
	std::vector<Vector> testVectors(numVectors, Vector(matrix.numberOfColumns()));

	testVectors[0] = tv;

	if (numVectors > 1) {
		Vector b(matrix.numberOfColumns(), 0.0);
#pragma omp parallel for
		for (count i = 1; i < numVectors; ++i) {
			for (count j = 0; j < matrix.numberOfColumns(); ++j) {
				testVectors[i][j] = 2 * Aux::Random::probability() - 1;
			}

			testVectors[i] = smoother.relax(matrix, b, testVectors[i], SETUP_TV_SWEEPS);
		}
	}

	return testVectors;
}

void MultiLevelSetup::addHighDegreeSeedNodes(const CSRMatrix &matrix, std::vector<index> &status) const {
	std::vector<count> deg(matrix.numberOfRows());
#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		deg[i] = matrix.nnzInRow(i) - 1;
	}

#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		double num = 0.0;
		double denom = 0.0;
		matrix.forNonZeroElementsInRow(i, [&](index j, double value){
			if (i != j) {
				num += std::abs(value) * (double) deg[j];
			} else {
				denom = std::abs(value);
			}
		});


		if ((double) deg[i] >= SETUP_AGGREGATION_DEGREE_THRESHOLD * (num / denom)) { // high degree node becomes seed
			status[i] = i;
		}
	}
}

void MultiLevelSetup::aggregateLooseNodes(const CSRMatrix &strongAdjMatrix, std::vector<index> &status, count &nc) const {
	std::vector<index> looseNodes;
	for (index i = 0; i < strongAdjMatrix.numberOfRows(); ++i) {
		double max = std::numeric_limits<double>::min();
		strongAdjMatrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (value > max) max = value;
		});

		if (std::abs(max) < 1e-9 || max == std::numeric_limits<double>::min()) {
			looseNodes.push_back(i);
		}
	}

	if (looseNodes.size() > 0) {
		status[looseNodes[0]] = looseNodes[0]; // mark first as seed
		for (index k = 1; k < looseNodes.size(); ++k) {
			status[looseNodes[k]] = looseNodes[0]; // first loose nodes becomes seed
		}

		nc -= looseNodes.size() - 1;
	}
}

void MultiLevelSetup::computeStrongAdjacencyMatrix(const CSRMatrix &matrix, CSRMatrix &strongAdjMatrix) const {
	std::vector<double> maxNeighbor(matrix.numberOfRows(), std::numeric_limits<double>::min());
#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
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
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		index cIdx = rowIdx[i];
		matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
				columnIdx[cIdx] = j;
				nonZeros[cIdx] = -value;
				++cIdx;
			}
		});
	}

	strongAdjMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), rowIdx, columnIdx, nonZeros, matrix.sorted());
}

void MultiLevelSetup::computeAffinityMatrix(const CSRMatrix &matrix, const std::vector<Vector> &tVs, CSRMatrix &affinityMatrix) const {
	assert(tVs.size() > 0);

	std::vector<index> rowIdx(matrix.numberOfRows()+1);
	std::vector<index> columnIdx(matrix.nnz());
	std::vector<double> nonZeros(matrix.nnz());

#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		rowIdx[i+1] = matrix.nnzInRow(i);
	}

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		rowIdx[i+1] += rowIdx[i];
	}

	std::vector<double> normSquared(matrix.numberOfRows(), 0.0);
#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		for (index k = 0; k < tVs.size(); ++k) {
			normSquared[i] += tVs[k][i] * tVs[k][i];
		}
	}

#pragma omp parallel for
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
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

	affinityMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), rowIdx, columnIdx, nonZeros, matrix.sorted());
}

void MultiLevelSetup::aggregationStage(const CSRMatrix &matrix, count &nc, const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, std::vector<Vector> &tVs, std::vector<index> &status) const {
	std::vector<std::vector<index>> bins(10);
	computeStrongNeighbors(affinityMatrix, status, bins);

	std::vector<double> diag(matrix.numberOfRows(), 0.0);
#pragma omp parallel for
	for (index i = 0 ; i < matrix.numberOfRows(); ++i) {
		diag[i] = matrix(i,i);
	}

	for (index k = bins.size(); k-- > 0;) { // iterate over undecided nodes with strong neighbors in decreasing order of strongest neighbor
		for (index i : bins[k]) {
			if (status[i] == UNDECIDED) { // node is still undecided
				index s = 0;
				if (findBestSeedEnergyCorrected(strongAdjMatrix, affinityMatrix, diag, tVs, status, i, s)) {
					status[s] = s; // s becomes seed
					status[i] = s; // i's seed is s
					nc--;

					for (index j = 0; j < tVs.size(); ++j) { // update test vectors
						tVs[j][i] = tVs[j][s];
					}
				}
			}
		}

		if (nc <= matrix.numberOfRows() * SETUP_COARSENING_WORK_GUARD / SETUP_CYCLE_INDEX) {
			break;
		}
	} // iterate over bins
}

void MultiLevelSetup::computeStrongNeighbors(const CSRMatrix &affinityMatrix, const std::vector<index> &status, std::vector<std::vector<index>> &bins) const {
	std::vector<bool> undecided(affinityMatrix.numberOfRows(), false);
	std::vector<double> maxNeighbor(affinityMatrix.numberOfRows(), std::numeric_limits<double>::min());
	double overallMax = 0.0;
	double overallMin = std::numeric_limits<double>::max();

	affinityMatrix.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) { // determine the highest affinity neighbor of each node
		if (status[i] == UNDECIDED && (status[j] == UNDECIDED || status[j] == j)) { // i is UNDECIDED and its neighbor j is also UNDECIDED or SEED
			if (value > maxNeighbor[i]) {
				maxNeighbor[i] = value;
			}
			undecided[i] = true;
		}
	});

	for (index i = 0; i < affinityMatrix.numberOfRows(); ++i) {
		if (maxNeighbor[i] > overallMax) {
			overallMax = maxNeighbor[i];
		}
		if (maxNeighbor[i] < overallMin) {
			overallMin = maxNeighbor[i];
		}
	}

	double h = fabs(overallMax - overallMin) < 1e-15? 1.0 : (double) bins.size() / (overallMax - overallMin);
	for (index i = 0; i < affinityMatrix.numberOfRows(); ++i) {
		if (undecided[i]) { // undecided nodes with strong neighbors
			index binIndex = (index) std::floor(h * (maxNeighbor[i] - overallMin));
			if (binIndex == bins.size()) { // last interval is closed on the right
				binIndex--;
			}

			assert(binIndex >= 0 && binIndex < bins.size());
			bins[binIndex].push_back(i);
		}
	}
}

bool MultiLevelSetup::findBestSeedEnergyCorrected(const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, const std::vector<double> &diag, const std::vector<Vector> &tVs, const std::vector<index> &status, const index u, index &s) const {
	bool foundSeed = false;
	std::vector<double> r(tVs.size(), 0.0);
	std::vector<double> q(tVs.size(), 0.0);
	std::vector<double> E(tVs.size(), 0.0);

	double d = diag[u];
	double d2 = 0.5 * diag[u];
	for (index k = 0; k < tVs.size(); ++k) {
		double rr = 0.0;
		double qq = 0.0;
		strongAdjMatrix.forNonZeroElementsInRow(u, [&](index v, double value) {
			rr += value * tVs[k][v];
			qq += value * 0.5 * tVs[k][v] * tVs[k][v];
		});

		r[k] = rr;
		q[k] = qq;
		double y = rr/d;
		E[k] = (d2*y - rr)*y + qq;
	}

	double maxNeighbor = -1.0;
	affinityMatrix.forNonZeroElementsInRow(u, [&](index v, double value) {
		if (status[v] == UNDECIDED || status[v] == v) {
			double maxMu = -1.0;
			bool smallRatio = true;
			for (index k = 0; k < tVs.size(); ++k) {
				double xv = tVs[k][v];
				double Ec = (d2*xv - r[k])*xv + q[k];
				double mu = Ec / (E[k] + 1e-15);

				if (mu > maxMu) {
					maxMu = mu;
				}
				if (maxMu > 2.5) {
					smallRatio = false;
					break;
				}
			}

			if (smallRatio && value > maxNeighbor) {
				maxNeighbor = value;
				s = v;
				foundSeed = true;
			}
		}
	});

	return foundSeed;
}

bool MultiLevelSetup::canCoarsen(const CSRMatrix &A) const {
	return A.numberOfRows() > MAX_DIRECT_SOLVE_SIZE;
}

bool MultiLevelSetup::isRelaxationFast(const CSRMatrix &A, index lvlIndex, Vector &tv) const {
	count nu = SETUP_RELAX_ACF_MIN_SWEEPS + 2 * (lvlIndex - 1);
	count tvNu = SETUP_TV_SWEEPS;
	count initial = 3;

	// create testVector in [-1,1]
	tv = Vector(A.numberOfRows());
	for (index i = 0; i < tv.getDimension(); ++i) {
		tv[i] = 2.0 * Aux::Random::probability() - 1.0;
	}

	Vector b(A.numberOfRows(), 0.0);
	Vector x = tv;
	x = smoother.relax(A, b, x, initial);
	tv = smoother.relax(A, b, x, tvNu - initial);
	Vector y = smoother.relax(A, b, tv, nu - tvNu);
	double relaxAcf = std::pow((y - y.mean()).length() / (x - x.mean()).length(), (double) 1.0 / (double) (nu - initial));
	return relaxAcf <= SETUP_MAX_COARSE_RELAX_ACF || !canCoarsen(A);
}

void MultiLevelSetup::galerkinOperator(const CSRMatrix &P, const CSRMatrix &A, const std::vector<index> &PColIndex, const std::vector<std::vector<index>> &PRowIndex, CSRMatrix &B) const {
	std::vector<CSRMatrix::Triple> triples;
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
			triples.push_back({i,j,value});
		});

		spa.increaseRow();
	}

	B = CSRMatrix(P.numberOfColumns(), P.numberOfColumns(), triples, true);
}


} /* namespace NetworKit */
