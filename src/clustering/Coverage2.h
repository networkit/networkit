/*
 * Coverage2.h
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef Coverage2_H_
#define Coverage2_H_

#include "QualityMeasure2.h"

namespace NetworKit {
	
/**
 * Coverage is the fraction of intra-cluster edges.
 */
class Coverage2: public NetworKit::QualityMeasure2 {
public:
	Coverage2();
	virtual ~Coverage2();

	virtual double getQuality(const Partition& zeta, const Graph& G);

};

} /* namespace NetworKit */
#endif /* Coverage2_H_ */
