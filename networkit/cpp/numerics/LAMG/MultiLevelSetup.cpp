/*
 * MultiLevelSetup.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#include "MultiLevelSetup.h"
#include "Level/LevelElimination.h"
#include "Level/EliminationStage.h"
#include "LAMGSettings.h"
#include "../../algebraic/LaplacianMatrix.h"
#include "../../io/LineFileReader.h"
#include <fstream>
#include <iostream>
#include <sstream>

#include "../../auxiliary/Enforce.h"
#include "../../auxiliary/Timer.h"
#include "../../algebraic/CSRMatrix.h"

#include <cstdio>
#include <set>

namespace NetworKit {

MultiLevelSetup::MultiLevelSetup(const Smoother &smoother) : smoother(smoother) {
}

void MultiLevelSetup::setup(const Graph &G, LevelHierarchy &hierarchy) const {
	setup(CSRMatrix::graphLaplacian(G), hierarchy);
}

void MultiLevelSetup::setup(const CSRMatrix &matrix, LevelHierarchy &hierarchy) const {
	CSRMatrix A = matrix;
	hierarchy.addFinestLevel(A);
	DEBUG("FINEST\t", matrix.numberOfRows(), "\t", matrix.nnz() - matrix.numberOfRows());

	bool doneCoarsening = false;
	count numTVs = TV_NUM;
	index level = 0;
	while (!doneCoarsening) {
		// ELIMINATION
		if (coarseningElimination(A, hierarchy)) {
			if (!canCoarsen(A)) doneCoarsening = true;
			level++;
			DEBUG(level, " ELIM\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
		}

		// AGGREGATION
		Vector tv;
		if (doneCoarsening || !isRelaxationFast(A, level, tv)) {
			doneCoarsening = true;
		} else {
			coarseningAggregation(A, hierarchy, tv, numTVs);
			level++;
			DEBUG(level, " AGG\t\t", A.numberOfRows(), "\t", A.nnz() / 2);
			if (numTVs < TV_MAX) {
				numTVs += TV_INC;
			}
		}

		if (!canCoarsen(A)) doneCoarsening = true;
	}
}

bool MultiLevelSetup::coarseningElimination(CSRMatrix &matrix, LevelHierarchy &hierarchy) const {
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
		for (index i = 0, fIndex = 0, cIndex = 0; i < matrix.numberOfRows(); ++i) {
			if (fNode[i] && fIndex < nf) {
				coarseIndex[i] = fIndex;
				fSet[fIndex++] = i;
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

		CSRMatrix Acc; // Schur complement
		CSRMatrix Acf; // Schur complement

		subMatrix(matrix, cSet, cSet, coarseIndex, Acc);
		subMatrix(matrix, cSet, fSet, coarseIndex, Acf);

		matrix = Acc + Acf * P;

		// check that sum of each row is zero (matrix is laplacian)
		for (index i = 0; i < matrix.numberOfRows(); ++i) {
			double sum = 0.0;
			matrix.forNonZeroElementsInRow(i, [&](index j, double value){
				sum += value;
			});

			sum -= matrix(i,i);
			matrix.setValue(i, i, -sum);
		}

		stageNum++;

		DEBUG("Elimination stage ", stageNum, ": total=", nc, " f=", nf, " c=", nc);
	}

	if (stageNum != 0) { // we have coarsened the matrix
		hierarchy.addEliminationLevel(matrix, coarseningStages);
		return true;
	}

	return false;
}

count MultiLevelSetup::lowDegreeSweep(const CSRMatrix &matrix, std::vector<bool> &fNode, index stage) const {
	fNode.resize(matrix.numberOfRows(), true); // first mark all nodes as f nodes
	count numFNodes = 0;
	int degreeOffset = stage != 0;

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		if ((int) matrix.nnzInRow(i) - degreeOffset <= SETUP_ELIMINATION_MAX_DEGREE && fNode[i]) { // node i has degree <= 4 and can be eliminated
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
	std::vector<std::pair<index, index>> positions;
	std::vector<double> values;
	q = Vector(fSet.size());
	for (index k = 0; k < fSet.size(); ++k) { // Afc
		matrix.forNonZeroElementsInRow(fSet[k], [&](index j, edgeweight w){
			if (fSet[k] == j) {
				q[k] = 1.0 / w;
			} else {
				positions.push_back(std::make_pair(k,coarseIndex[j]));
				values.push_back(w);
			}
		});
	}

	for (index i = 0; i < values.size(); ++i) { // * Aff^-1
		values[i] *= -q[positions[i].first];
	}

	P = CSRMatrix(fSet.size(), coarseIndex.size() - fSet.size(), positions, values);
}

void MultiLevelSetup::subMatrix(const CSRMatrix &matrix, const std::vector<index> &rows, const std::vector<index> &columns, const std::vector<index> &coarseIndex, CSRMatrix &result) const {
	std::vector<std::pair<index, index>> positions;
	std::vector<double> values;

	for (index k = 0; k < rows.size(); ++k) {
		matrix.forNonZeroElementsInRow(rows[k], [&](index j, edgeweight value) {
			if (coarseIndex[j] < columns.size() && columns[coarseIndex[j]] == j) { // check if neighbor is in columns
				positions.push_back(std::make_pair(k,coarseIndex[j]));
				values.push_back(value);
			}
		});
	}

	result = CSRMatrix(rows.size(), columns.size(), positions, values);
}

void MultiLevelSetup::coarseningAggregation(CSRMatrix &matrix, LevelHierarchy &hierarchy, Vector &tv, count numTVVectors) const {
	Vector B(SETUP_MAX_AGGREGATION_STAGES, std::numeric_limits<double>::max());
	std::vector<std::vector<int64_t>> S(SETUP_MAX_AGGREGATION_STAGES, std::vector<int64_t>(matrix.numberOfRows(), std::numeric_limits<int64_t>::max()));
	std::vector<int64_t> status(matrix.numberOfRows(), UNDECIDED);
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
		alpha <= maxCoarseningRatio? B[stage] = 1-alpha : B[stage] = 1+alpha;

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
	std::vector<std::pair<index, index>> pPositions(matrix.numberOfRows());
	std::vector<double> pValues(matrix.numberOfRows());
	std::vector<index> PColIndex(matrix.numberOfRows());
	std::vector<std::vector<index>> PRowIndex(nc[bestAggregate]);
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		pPositions[i] = std::make_pair(i, status[i]);
		pValues[i] = 1;
		PColIndex[i] = status[i];
		PRowIndex[status[i]].push_back(i);
	}

	CSRMatrix P(matrix.numberOfRows(), nc[bestAggregate], pPositions, pValues);

	// create coarsened laplacian
	galerkinOperator(P, matrix, PColIndex, PRowIndex, matrix);

	// check that sum of each row is zero (matrix is laplacian)
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		edgeweight sum = 0.0;
		matrix.forNonZeroElementsInRow(i, [&](index j, edgeweight value){
			sum += value;
		});

		sum -= matrix(i,i);
		matrix.setValue(i, i, -sum);
	}
	hierarchy.addAggregationLevel(matrix, P);
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

void MultiLevelSetup::addHighDegreeSeedNodes(const CSRMatrix &matrix, std::vector<int64_t> &status) const {
	std::vector<count> deg(matrix.numberOfRows());
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

void MultiLevelSetup::aggregateLooseNodes(const CSRMatrix &strongAdjMatrix, std::vector<int64_t> &status, count &nc) const {
	std::vector<index> looseNodes;
	for (index i = 0; i < strongAdjMatrix.numberOfRows(); ++i) {
		double max = std::numeric_limits<double>::min();
		strongAdjMatrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (value > max) max = value;
		});

		if (std::abs(max) < 1e-14 || max == std::numeric_limits<double>::min()) {
			looseNodes.push_back(i);
		}
	}

	if (looseNodes.size() > 0) {
		status[looseNodes[0]] = looseNodes[0]; // mark first as seed
		for (index k = 1; k < looseNodes.size(); ++k) {
			status[looseNodes[k]] = looseNodes[0]; // first loose nodes becomes seed
		}

		nc -= looseNodes.size() + 1;
	}
}

void MultiLevelSetup::computeStrongAdjacencyMatrix(const CSRMatrix &matrix, CSRMatrix &strongAdjMatrix) const {
	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;

	std::vector<double> maxNeighbor(matrix.numberOfRows(), std::numeric_limits<double>::min());
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (i != j) {
				if (-value > maxNeighbor[i]) {
					maxNeighbor[i] = -value;
				}
			}
		});
	}

	matrix.forNonZeroElementsInRowOrder([&](index i, index j, double value) {
		if (i != j && std::abs(value) >= 0.1 * std::min(maxNeighbor[i], maxNeighbor[j])) {
			positions.push_back(std::make_pair(i,j));
			values.push_back(-value);
		}
	});

	strongAdjMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), positions, values);
}

void MultiLevelSetup::computeAffinityMatrix(const CSRMatrix &matrix, const std::vector<Vector> &tVs, CSRMatrix &affinityMatrix) const {
	assert(tVs.size() > 0);

	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;

	std::vector<double> normSquared(matrix.numberOfRows(), 0.0);
	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		for (index k = 0; k < tVs.size(); ++k) {
			normSquared[i] += tVs[k][i] * tVs[k][i];
		}
	}

	for (index i = 0; i < matrix.numberOfRows(); ++i) {
		double nir = 1.0 / normSquared[i];
		matrix.forNonZeroElementsInRow(i, [&](index j, double value) {
			if (j >= i) {
				double ij = 0.0;
				for (index k = 0; k < tVs.size(); ++k) {
					ij += tVs[k][i] * tVs[k][j];
				}

				double value = (ij * ij) * nir / normSquared[j];
				positions.push_back(std::make_pair(i,j));
				values.push_back(value);

				if (j > i) { // symmetric
					positions.push_back(std::make_pair(j,i));
					values.push_back(value);
				}
			}
		});
	}

	affinityMatrix = CSRMatrix(matrix.numberOfRows(), matrix.numberOfColumns(), positions, values);
}

void MultiLevelSetup::aggregationStage(const CSRMatrix &matrix, count &nc, const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, std::vector<Vector> &tVs, std::vector<int64_t> &status) const {
	std::vector<std::vector<index>> bins(10);
	computeStrongNeighbors(affinityMatrix, status, bins);

	std::vector<double> diag(matrix.numberOfRows(), 0.0);
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

void MultiLevelSetup::computeStrongNeighbors(const CSRMatrix &affinityMatrix, const std::vector<int64_t> &status, std::vector<std::vector<index>> &bins) const {
	std::vector<bool> undecided(affinityMatrix.numberOfRows(), false);
	std::vector<double> maxNeighbor(affinityMatrix.numberOfRows(), std::numeric_limits<double>::min());
	double overallMax = 0.0;
	double overallMin = std::numeric_limits<double>::max();

	affinityMatrix.forNonZeroElementsInRowOrder([&](index i, index j, double value) { // determine the highest affinity neighbor of each node
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

bool MultiLevelSetup::findBestSeed(const CSRMatrix &affinityMatrix, const std::vector<index> &strongNeighborsOfU, const std::vector<int64_t> &status, const index u, index &s) const {
	double maxAffinity = std::numeric_limits<double>::min();
	for (index i = 0; i < strongNeighborsOfU.size(); ++i) {
		index v = strongNeighborsOfU[i];
		if (status[v] < 0 || (index) status[v] == v) { // neighbor is seed or undecided
			if (affinityMatrix(u, v) > maxAffinity) {
				s = v;
				maxAffinity = affinityMatrix(u, v);
			}
		}
	}

	return maxAffinity != std::numeric_limits<double>::min(); // we have found a strong neighbor which is seed or undecided
}

bool MultiLevelSetup::findBestSeedEnergyCorrected(const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, const std::vector<double> &diag, const std::vector<Vector> &tVs, const std::vector<int64_t> &status, const index u, index &s) const {
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

	return relaxAcf > SETUP_MAX_COARSE_RELAX_ACF;
}

void MultiLevelSetup::galerkinOperator(const CSRMatrix &P, const CSRMatrix &A, const std::vector<index> &PColIndex, const std::vector<std::vector<index>> &PRowIndex, CSRMatrix &B) const {
	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;
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
			positions.push_back(std::make_pair(i,j));
			values.push_back(value);
		});

		spa.increaseRow();
	}

	B = CSRMatrix(P.numberOfColumns(), P.numberOfColumns(), positions, values);
}


} /* namespace NetworKit */
