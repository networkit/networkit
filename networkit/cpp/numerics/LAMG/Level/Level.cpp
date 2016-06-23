/*
 * Level.cpp
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#include "Level.h"
#include "../LAMGSettings.h"

namespace NetworKit {

Level::Level(LevelType type) : type(type) {
}

Level::Level(LevelType type, const CSRMatrix &A) : type(type), A(A) {
}

const CSRMatrix& Level::getLaplacian() const {
	return A;
}

count Level::getNumberOfNodes() const {
	return A.numberOfRows();
}

void Level::restrict(const Vector &bf, Vector &bc) const {
}

void Level::restrict(const Vector &bf, Vector &bc, std::vector<Vector> &bStages) const {
}

void Level::interpolate(const Vector &xc, Vector &xf) const {
}

void Level::interpolate(const Vector &xc, Vector &xf, const std::vector<Vector> &bStages) const {
}



} /* namespace NetworKit */
