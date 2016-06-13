/*
 * LevelHierarchy.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#include "LevelHierarchy.h"
#include "LAMGSettings.h"

namespace NetworKit {

LevelHierarchy::LevelHierarchy() {
}

void LevelHierarchy::addFinestLevel(const CSRMatrix &A) {
	finestLevel = LevelFinest(A);
}

void LevelHierarchy::addEliminationLevel(const CSRMatrix &A, const std::vector<EliminationStage> &coarseningStages) {
	levelType.push_back(ELIMINATION);
	levelIndex.push_back(eliminationLevels.size());
	eliminationLevels.push_back(LevelElimination(A, coarseningStages));
}

void LevelHierarchy::addAggregationLevel(const CSRMatrix &A, const CSRMatrix &P, const CSRMatrix &R) {
	levelType.push_back(AGGREGATION);
	levelIndex.push_back(aggregationLevels.size());
	aggregationLevels.push_back(LevelAggregation(A, P, R));
}

void LevelHierarchy::setLastAsCoarsest() {
	CSRMatrix A = this->at(size()-1).getLaplacian();
	count n = A.numberOfRows() + 1;
	std::vector<double> entries(n*n, 0.0);
	A.parallelForNonZeroElementsInRowOrder([&](index i, index j, double value) {
		entries[i * n + j] = value;
	});

	for (index i = 0; i < n-1; ++i) {
		entries[i * n + n-1] = 1;
		entries[(n-1)*n + i] = 1;
	}

	coarseLUMatrix = DenseMatrix(n, n, entries);
	DenseMatrix::LUDecomposition(coarseLUMatrix);
}

DenseMatrix& LevelHierarchy::getCoarseMatrix() {
	return coarseLUMatrix;
}

count LevelHierarchy::size() const {
	return levelType.size() + 1; // elimination + aggregation levels + finestLevel
}

LevelType LevelHierarchy::getType(index levelIdx) const {
	assert(levelIdx >= 0 && levelIdx < this->size());

	if (levelIdx == 0) {
		return FINEST;
	} else {
		return levelType[levelIdx-1];
	}
}

Level& LevelHierarchy::at(index levelIdx) {
	assert(levelIdx >= 0 && levelIdx < this->size());

	if (levelIdx == 0) { // finest level
		return finestLevel;
	} else {
		levelIdx--;
	}

	if (levelType[levelIdx] == ELIMINATION) {
		return eliminationLevels[levelIndex[levelIdx]];
	} else {
		return aggregationLevels[levelIndex[levelIdx]];
	}
}

double LevelHierarchy::cycleIndex(index levelIdx) {
	double gamma = 1.0;
	if (getType(levelIdx+1) != ELIMINATION) {
		count finestNumEdges = finestLevel.getLaplacian().nnz();
		count numFineEdges = this->at(levelIdx).getLaplacian().nnz();

		if (numFineEdges > 0.1 * finestNumEdges) {
			gamma = SETUP_CYCLE_INDEX;
		} else {
			count numCoarseEdges = this->at(levelIdx+1).getLaplacian().nnz();
			gamma = std::max(1.0, std::min(1.5, SETUP_COARSENING_WORK_GUARD / ((double) numCoarseEdges / (double) numFineEdges)));
		}
	}

	return gamma;
}

} /* namespace NetworKit */
