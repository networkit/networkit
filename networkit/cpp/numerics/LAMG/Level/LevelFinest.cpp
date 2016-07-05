/*
 * LevelFinest.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#include "LevelFinest.h"

namespace NetworKit {

LevelFinest::LevelFinest() : Level(FINEST) {
}

LevelFinest::LevelFinest(const CSRMatrix &A) : Level(LevelType::FINEST, A) {
	// set number of testVectors to TV_NUM??
}

void LevelFinest::coarseType(const Vector &xf, Vector &xc) const {
	// do nothing!
}

void LevelFinest::restrict(const Vector &bf, Vector &bc) const {
	// do nothing!
}

void LevelFinest::interpolate(const Vector &xc, Vector &xf) const {
	// do nothing!
}

} /* namespace NetworKit */
