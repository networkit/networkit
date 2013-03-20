/*
 * Coverage.h
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef COVERAGE_H_
#define COVERAGE_H_

#include "QualityMeasure.h"

namespace EnsembleClustering {

class Coverage: public EnsembleClustering::QualityMeasure {
public:
	Coverage();
	virtual ~Coverage();

	virtual double getQuality(const Clustering& zeta, Graph& G);

};

} /* namespace EnsembleClustering */
#endif /* COVERAGE_H_ */
