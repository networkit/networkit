/*
 * MultigridHierarchy.cpp
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultigridHierarchy.h"

namespace NetworKit {

MultigridHierarchy::MultigridHierarchy(const Matrix &fineMatrix) {
	laplacianMatrices.push_back(fineMatrix);
}

void MultigridHierarchy::addLevel(const Matrix &coarseMatrix, const Matrix &interpolationMatrix) {
	laplacianMatrices.push_back(coarseMatrix);
	interpolationMatrices.push_back(interpolationMatrix);
}

const Matrix& MultigridHierarchy::getLaplacian(const index level) {
	assert(level < getNumLevels());
	return laplacianMatrices[level];
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

float MultigridHierarchy::getNumMultigridCycles(const index level) const {
	return 0.0;
}

Vector MultigridHierarchy::createInitialCoarseResult(const index level, const Vector &xFine) const {
	return Vector(laplacianMatrices[level].numberOfRows());
}

Vector MultigridHierarchy::restriction(const index level, const Vector &xFine, const Vector &bFine) const {
	assert(level < getNumLevels());
	Vector residual = bFine - laplacianMatrices[level] * xFine;
	return Matrix::mTvMultiply(interpolationMatrices[level], residual);
}

Vector MultigridHierarchy::prolongation(const index level, const Vector &xCoarse, const Vector &xFine, const Vector &bFine) const {
	assert(level < getNumLevels());
	return xFine + interpolationMatrices[level] * xCoarse;
}

count MultigridHierarchy::getNumLevels() const {
	return laplacianMatrices.size();
}

} /* namespace NetworKit */
