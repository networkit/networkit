/*
 * LevelAggregation.cpp
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#include "LevelAggregation.h"

namespace NetworKit {

LevelAggregation::LevelAggregation(const CSRMatrix &A, const CSRMatrix &P, const CSRMatrix &R) : Level(LevelType::AGGREGATION, A), P(P), R(R) {
}

void LevelAggregation::coarseType(const Vector &xf, Vector &xc) const {
	xc = Vector(P.numberOfColumns(), 0.0);
}

void LevelAggregation::restrict(const Vector &bf, Vector &bc) const {
	bc = R * bf;
}

void LevelAggregation::interpolate(const Vector &xc, Vector &xf) const {
	xf = P * xc;
}

} /* namespace NetworKit */
