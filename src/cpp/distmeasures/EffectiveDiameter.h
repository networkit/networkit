/*
 * EffectiveDiameter.h
 *
 *  Created on: Jun 14, 2014
 *      Author: Marc Nemes
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

	EffectiveDiameter();

	virtual ~EffectiveDiameter();

	/**
	 * determines the effective diameter of a given graph
	 */
	int effectiveDiameter(const Graph& G);

};
}/* namespace NetworKit */


















#endif /* EFFECTIVEDIAMETER_H_ */
