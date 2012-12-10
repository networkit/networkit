/*
 * QualityMeasure.h
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#ifndef QUALITYMEASURE_H_
#define QUALITYMEASURE_H_

#include "Clustering.h"

namespace EnsembleClustering {

/**
 * Abstract base class for all clustering quality measures.
 */
class QualityMeasure {

protected:

	Graph* G;

public:

	QualityMeasure(Graph& G);

	virtual ~QualityMeasure();

	virtual double getQuality(Clustering& zeta) =0;
};

} /* namespace EnsembleClustering */
#endif /* QUALITYMEASURE_H_ */
