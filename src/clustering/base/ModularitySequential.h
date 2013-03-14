/*
 * ModularitySequential.h
 *
 *  Created on: 14.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MODULARITYSEQUENTIAL_H_
#define MODULARITYSEQUENTIAL_H_

#include <unordered_map>
#include <cmath>
#include <stdexcept>

#include "QualityMeasure.h"
#include "CoverageSequential.h"

#include "../../base/IndexMap.h"
#include "../../graph/NodeMap.h"

namespace EnsembleClustering {

class ModularitySequential: public EnsembleClustering::QualityMeasure {
public:
	ModularitySequential();
	virtual ~ModularitySequential();

	virtual double getQuality(const Clustering& zeta, Graph& G);

};

} /* namespace EnsembleClustering */
#endif /* MODULARITYSEQUENTIAL_H_ */
