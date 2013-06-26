/*
 * Surprise.h
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#ifndef SURPRISE_H_
#define SURPRISE_H_

#include "QualityMeasure.h"

#include <cmath>
#include <algorithm>

#include "../auxiliary/MissingMath.h"

namespace NetworKit {

class Surprise: public NetworKit::QualityMeasure {

public:

	Surprise();

	virtual ~Surprise();

	virtual double getQuality(const Clustering& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif /* SURPRISE_H_ */
