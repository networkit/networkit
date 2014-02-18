/*
 * ModularitySequential.h
 *
 *  Created on: 14.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MODULARITYSEQUENTIAL_H_
#define MODULARITYSEQUENTIAL_H_

#include "QualityMeasure.h"

namespace NetworKit {

class ModularitySequential: public NetworKit::QualityMeasure {
public:
	ModularitySequential();
	virtual ~ModularitySequential();

	virtual double getQuality(const Partition& zeta, const Graph& G);

};

} /* namespace NetworKit */
#endif /* MODULARITYSEQUENTIAL_H_ */
