/*
 * MultigridHierarchy.cpp
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "MultigridHierarchy.h"

namespace NetworKit {

void MultigridHierarchy::addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix) {
	laplacians.push_back(laplacian);
	interpolationMatrices.push_back(interpolationMatrix);
}

const Matrix& MultigridHierarchy::getLaplacian(const index level) {
	assert(level < getNumLevels());
	return laplacians[level];
}

const Matrix& MultigridHierarchy::getInterpolationMatrix(const index level) {
	assert(level < getNumLevels());
	return interpolationMatrices[level];
}

count MultigridHierarchy::getNumLevels() const {
	return laplacians.size();
}

} /* namespace NetworKit */
