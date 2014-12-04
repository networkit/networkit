/*
 * LAMGHierarchy.cpp
 *
 *  Created on: 14.11.2014
 *      Author: Michael
 */

#include "LAMGHierarchy.h"

namespace NetworKit {

LAMGHierarchy::LAMGHierarchy(const Matrix &fineMatrix) : MultigridHierarchy(fineMatrix) {
}

void LAMGHierarchy::addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, Type type) {
	addLevel(laplacian, interpolationMatrix);
	this->type.push_back(type);
}

LAMGHierarchy::Type LAMGHierarchy::getType(const index level) const {
	assert(level < getNumLevels());
	return type[level];
}

count LAMGHierarchy::nnzOfMatrix(const index level) const {
	assert(level < getNumLevels());
	return coarseMatrices[level].nnz();
}

count LAMGHierarchy::getNumMultigridCycles(const index level) const {
	assert(level < getNumLevels() - 1);
	count cycles = 1;
	if (type[level + 1] != ELIMINATION) {
		count nnzFinest = nnzOfMatrix(0);
		count nnzCurrent = nnzOfMatrix(level);
		if (nnzCurrent > 0.1 * nnzFinest) {
			if (Aux::Random::real(1) <= 0.5) cycles++;
		} else {
			count nnzNext = nnzOfMatrix(level + 1);
			double ratio = std::min(2.0, 0.7 * (double) nnzNext / (double) nnzCurrent);
			if (ratio > 1) {
				if (Aux::Random::real(1) <= 2 - ratio) cycles++;
			}
		}
	}

	return cycles;
}

count LAMGHierarchy::getNumPreSmooth(const index level) const {
	assert(level < getNumLevels() - 1);
	if (type[level + 1] == ELIMINATION) {
		return 0;
	} else {
		return 1;
	}
}

count LAMGHierarchy::getNumPostSmooth(const index level) const {
	assert(level < getNumLevels() - 1);
	if (type[level + 1] == ELIMINATION) {
		return 0;
	} else {
		return 2;
	}
}

Vector LAMGHierarchy::restriction(const index level, const Vector &currentApproximation, const Vector &rhs) const {
	assert(level < getNumLevels());
	Vector restrictedVector;

	switch(type[level]) {
		case ELIMINATION:
			restrictedVector = Matrix::mTvMultiply(interpolationMatrices[level], rhs);
			break;
		case AGGREGATION:
			restrictedVector = Matrix::mTvMultiply(interpolationMatrices[level], rhs - coarseMatrices[level] * currentApproximation);
			break;
		default:
			break;
	}

	return restrictedVector;
}

Vector LAMGHierarchy::prolongation(const index level, const Vector &coarseResult, const Vector &fineResult) const {
	assert(level < getNumLevels());

	Vector prolongatedVector;
	switch(type[level]) {
		case ELIMINATION:
			prolongatedVector = interpolationMatrices[level] * coarseResult;
			break;
		case AGGREGATION:
			prolongatedVector = fineResult + interpolationMatrices[level] * coarseResult;
			break;
		default:
			break;
	}

	return prolongatedVector;
}

} /* namespace NetworKit */
