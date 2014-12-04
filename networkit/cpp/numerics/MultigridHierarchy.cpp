/*
 * MultigridHierarchy.cpp
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultigridHierarchy.h"

namespace NetworKit {

MultigridHierarchy::MultigridHierarchy(const Matrix &fineMatrix) {
	coarseMatrices.push_back(fineMatrix);
}

void MultigridHierarchy::addLevel(const Matrix &coarseMatrix, const Matrix &interpolationMatrix) {
	coarseMatrices.push_back(coarseMatrix);
	interpolationMatrices.push_back(interpolationMatrix);
}

const Matrix& MultigridHierarchy::getLaplacian(const index level) {
	assert(level < getNumLevels());
	return coarseMatrices[level];
}

const Matrix& MultigridHierarchy::getInterpolationMatrix(const index level) {
	assert(level < getNumLevels());
	return interpolationMatrices[level];
}

count MultigridHierarchy::getNumPreSmooth(const index level) const {
	return 3;
}

count MultigridHierarchy::getNumPostSmooth(const index level) const {
	return 3;
}

count MultigridHierarchy::getNumMultigridCycles(const index level) const {
	return 0;
}

Vector MultigridHierarchy::restriction(const index level, const Vector &currentApproximation, const Vector &rhs) const {
	assert(level < getNumLevels());
	Vector residual = rhs - coarseMatrices[level] * currentApproximation;
	return Matrix::mTvMultiply(interpolationMatrices[level], residual);
}

Vector MultigridHierarchy::prolongation(const index level, const Vector &coarseApproximation, const Vector &fineApproximation) const {
	assert(level < getNumLevels());
	return fineApproximation + interpolationMatrices[level] * coarseApproximation;
}

count MultigridHierarchy::getNumLevels() const {
	return coarseMatrices.size();
}

} /* namespace NetworKit */
