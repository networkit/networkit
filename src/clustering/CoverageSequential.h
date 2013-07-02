/*
 * CoverageSequential.h
 *
 *  Created on: 14.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef COVERAGESEQUENTIAL_H_
#define COVERAGESEQUENTIAL_H_

#include "QualityMeasure.h"

namespace NetworKit {

class CoverageSequential: public NetworKit::QualityMeasure {
public:
	CoverageSequential();
	virtual ~CoverageSequential();

	virtual double getQuality(const Clustering& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif /* COVERAGESEQUENTIAL_H_ */
