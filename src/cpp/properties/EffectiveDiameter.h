/*
 * Eccentricity.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef EFFECTIVEDIAMETER_H_
#define EFFECTIVEDIAMETER_H_

#include "../graph/Graph.h"

#include <vector>
#include <set>
#include <math.h>
#include <iterator>

namespace NetworKit {

class EffectiveDiameter {

public:
	/**
	 * determines the effective diameter of a given graph
	 */
	static int effectiveDiameter(const Graph& G);

};

} /* namespace NetworKit */

#endif /* ECCENTRICITY_H_ */
