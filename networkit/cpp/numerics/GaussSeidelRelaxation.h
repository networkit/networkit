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

class GaussSeidelRelaxation : public Smoother {

private:
	double tolerance;

public:
	GaussSeidelRelaxation(double tolerance=1e-6);

	Vector relax(const Matrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations = std::numeric_limits<count>::max()) const;

};

} /* namespace NetworKit */

#endif /* GAUSSSEIDELRELAXATION_H_ */
