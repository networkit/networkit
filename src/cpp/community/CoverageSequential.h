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

/**
 * @deprecated DEPRECATED Sequential implementation for testing - use parallel implementation NetworKit::Coverage
 * 
 * Coverage is the fraction of intra-cluster edges.
 * 
 */
class CoverageSequential: public NetworKit::QualityMeasure {
public:
	/** Default destructor */
	virtual ~CoverageSequential();

	virtual double getQuality(const Partition& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif /* COVERAGESEQUENTIAL_H_ */
