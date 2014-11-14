/*
 * LAMGHierarchy.cpp
 *
 *  Created on: 14.11.2014
 *      Author: Michael
 */

#include "LAMGHierarchy.h"

namespace NetworKit {

void LAMGHierarchy::addLevel(const Matrix &laplacian, const Matrix &interpolationMatrix, Type type) {
	addLevel(laplacian, interpolationMatrix);
	this->type.push_back(type);
}

LAMGHierarchy::Type LAMGHierarchy::getType(const index level) const {
	assert(level < getNumLevels());
	return type[level];
}

} /* namespace NetworKit */
