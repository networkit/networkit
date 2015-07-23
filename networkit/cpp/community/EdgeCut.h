/*
 * EdgeCut.h
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#ifndef EDGECUT_H_
#define EDGECUT_H_

#include "QualityMeasure.h"

namespace NetworKit {

/**
 * @ingroup community
 */
class EdgeCut: public NetworKit::QualityMeasure {
public:
	virtual double getQuality(const Partition& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif /* EDGECUT_H_ */
