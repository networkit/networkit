/*
 * LevelFinest.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef LEVELFINEST_H_
#define LEVELFINEST_H_

#include "Level.h"

namespace NetworKit {

/**
 * @ingroup numerics
 */
class LevelFinest : public Level {
public:
	LevelFinest();
	LevelFinest(const CSRMatrix &A);

	void coarseType(const Vector &xf, Vector &xc) const;
	void restrict(const Vector &bf, Vector &bc);
	void interpolate(const Vector &xc, Vector &xf) const;
};

} /* namespace NetworKit */

#endif /* LEVELFINEST_H_ */
