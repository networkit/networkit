/*
 * QualityMeasure2.h
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef QualityMeasure2_H_
#define QualityMeasure2_H_

#include "../structures/Partition.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Abstract base class for all clustering quality measures.
 */
class QualityMeasure2 {


public:

	QualityMeasure2();

	virtual ~QualityMeasure2();

	virtual double getQuality(const Partition& zeta, const Graph& G) =0;
};

} /* namespace NetworKit */
#endif /* QualityMeasure2_H_ */
