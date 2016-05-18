/*
 * GaussSeidelRelaxation.h
 *
 *  Created on: 27.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef GAUSSSEIDELRELAXATION_H_
#define GAUSSSEIDELRELAXATION_H_

#include "Smoother.h"

namespace NetworKit {

/**
 * @ingroup numerics
 * Implementation of the Gauss-Seidel smoother.
 */
class GaussSeidelRelaxation : public Smoother {

private:
	double tolerance;

public:
	GaussSeidelRelaxation(double tolerance=1e-15);

	Vector relax(const CSRMatrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations = std::numeric_limits<count>::max()) const;
	Vector relax(const CSRMatrix &A, const Vector &b, const count maxIterations = std::numeric_limits<count>::max()) const;

};

} /* namespace NetworKit */

#endif /* GAUSSSEIDELRELAXATION_H_ */
