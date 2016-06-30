/*
 * LevelAggregation.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef LEVELAGGREGATION_H_
#define LEVELAGGREGATION_H_

#include "Level.h"

namespace NetworKit {

/**
 * @ingroup numerics
 */
class LevelAggregation : public Level {
private:
	CSRMatrix P; // interpolation matrix (n x nc)
	CSRMatrix R; // restriction matrix (nc x n)

public:
	LevelAggregation(const CSRMatrix &A, const CSRMatrix &P, const CSRMatrix &R);

	void coarseType(const Vector &xf, Vector &xc) const;

	void restrict(const Vector &bf, Vector &bc) const;

	void interpolate(const Vector &xc, Vector &xf) const;

};

} /* namespace NetworKit */

#endif /* LEVELAGGREGATION_H_ */
