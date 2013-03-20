/*
 * CoverageSequential.h
 *
 *  Created on: 14.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef COVERAGESEQUENTIAL_H_
#define COVERAGESEQUENTIAL_H_

#include "QualityMeasure.h"

namespace EnsembleClustering {

class CoverageSequential: public EnsembleClustering::QualityMeasure {
public:
	CoverageSequential();
	virtual ~CoverageSequential();

	virtual double getQuality(const Clustering& zeta, Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* COVERAGESEQUENTIAL_H_ */
