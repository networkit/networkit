/*
 * LevelHierarchy.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef LEVELHIERARCHY_H_
#define LEVELHIERARCHY_H_

#include "Level/Level.h"
#include "Level/LevelFinest.h"
#include "Level/LevelElimination.h"
#include "Level/LevelAggregation.h"
#include "LAMGSettings.h"
#include "../../algebraic/DenseMatrix.h"

namespace NetworKit {

/**
 * @ingroup numerics
 */
template<class Matrix>
class LevelHierarchy {
private:
	std::vector<LevelType> levelType;
	std::vector<index> levelIndex;
	std::vector<LevelElimination<Matrix>> eliminationLevels;
	std::vector<LevelAggregation<Matrix>> aggregationLevels;
	LevelFinest<Matrix> finestLevel;
	DenseMatrix coarseLUMatrix;

	void createCoarseMatrix();

public:
	LevelHierarchy() = default;

	void addFinestLevel(const Matrix& A);
	void addEliminationLevel(const Matrix& A, const std::vector<EliminationStage<Matrix>>& coarseningStages);
	void addAggregationLevel(const Matrix& A, const Matrix& P, const Matrix& R);
	void setLastAsCoarsest();
	inline DenseMatrix& getCoarseMatrix() {
		return coarseLUMatrix;
	}

	inline count size() const {
		return levelType.size() + 1; // elimination + aggregation levels + finestLevel
	}

	LevelType getType(index levelIdx) const;
	Level<Matrix>& at(index levelIdx);
	double cycleIndex(index levelIdx);
};

template<class Matrix>
void LevelHierarchy<Matrix>::addFinestLevel(const Matrix& A) {
	finestLevel = LevelFinest<Matrix>(A);
}

template<class Matrix>
void LevelHierarchy<Matrix>::addEliminationLevel(const Matrix& A, const std::vector<EliminationStage<Matrix>>& coarseningStages) {
	levelType.push_back(ELIMINATION);
	levelIndex.push_back(eliminationLevels.size());
	eliminationLevels.push_back(LevelElimination<Matrix>(A, coarseningStages));
}

template <class Matrix>
void LevelHierarchy<Matrix>::addAggregationLevel(const Matrix& A, const Matrix& P, const Matrix& R) {
	levelType.push_back(AGGREGATION);
	levelIndex.push_back(aggregationLevels.size());
	aggregationLevels.push_back(LevelAggregation<Matrix>(A, P, R));
}

template<class Matrix>
void LevelHierarchy<Matrix>::setLastAsCoarsest() {
	const Matrix& A = this->at(size()-1).getLaplacian();
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


template<class Matrix>
LevelType LevelHierarchy<Matrix>::getType(index levelIdx) const {
	assert(levelIdx >= 0 && levelIdx < this->size());

	if (levelIdx == 0) {
		return FINEST;
	} else {
		return levelType[levelIdx-1];
	}
}

template<class Matrix>
Level<Matrix>& LevelHierarchy<Matrix>::at(index levelIdx) {
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

template<class Matrix>
double LevelHierarchy<Matrix>::cycleIndex(index levelIdx) {
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

#endif /* LEVELHIERARCHY_H_ */
