/*
 * NearlyLinearLaplacianSmoother.h
 *
 *  Created on: Jun 6, 2015
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_NUMERICS_NEARLYLINEARLAPLACIANSMOOTHER_H_
#define NETWORKIT_CPP_NUMERICS_NEARLYLINEARLAPLACIANSMOOTHER_H_

#include "Smoother.h"

namespace NetworKit {

class NearlyLinearLaplacianSmoother : public Smoother {
private:
	double tolerance;
public:
	NearlyLinearLaplacianSmoother(double tolerance = 1e-8);
	~NearlyLinearLaplacianSmoother() = default;


	Vector relax(const CSRMatrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations = std::numeric_limits<count>::max()) const;
	Vector relax(const CSRMatrix &A, const Vector &b, const count maxIterations = std::numeric_limits<count>::max()) const;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_NUMERICS_NEARLYLINEARLAPLACIANSMOOTHER_H_ */
