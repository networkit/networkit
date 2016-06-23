/*
 * Smoother.h
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef SMOOTHER_H_
#define SMOOTHER_H_

#include "../algebraic/Matrix.h"
#include "../algebraic/CSRMatrix.h"
#include "../algebraic/Vector.h"

#include <limits>

namespace NetworKit {

/**
 * @ingroup numerics
 * Abstract base class of a smoother.
 */
class Smoother {
public:
	Smoother() {}
	virtual ~Smoother(){}

	virtual Vector relax(const CSRMatrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations = std::numeric_limits<count>::max()) const = 0;
	virtual Vector relax(const CSRMatrix &A, const Vector &b, const count maxIterations = std::numeric_limits<count>::max()) const = 0;
};

} /* namespace NetworKit */

#endif /* SMOOTHER_H_ */
