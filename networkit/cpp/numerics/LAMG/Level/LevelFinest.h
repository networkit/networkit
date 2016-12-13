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
template<class Matrix>
class LevelFinest : public Level<Matrix> {
public:
	LevelFinest() : Level<Matrix>(LevelType::FINEST) {}
	LevelFinest(const Matrix& A) : Level<Matrix>(LevelType::FINEST, A) {}
};

} /* namespace NetworKit */

#endif /* LEVELFINEST_H_ */
