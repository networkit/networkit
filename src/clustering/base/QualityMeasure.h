/*
 * QualityMeasure.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef QUALITYMEASURE_H_
#define QUALITYMEASURE_H_

#include "Clustering.h"
#include "../../graph/Graph.h"

namespace EnsembleClustering {

/**
 * Abstract base class for all clustering quality measures.
 */
class QualityMeasure {


public:

	QualityMeasure();

	virtual ~QualityMeasure();

	virtual double getQuality(const Clustering& zeta, Graph& G) =0;
};

} /* namespace EnsembleClustering */
#endif /* QUALITYMEASURE_H_ */
