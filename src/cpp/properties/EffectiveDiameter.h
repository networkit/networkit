/*
 * EffectiveDiameter.h
 *
 *  Created on: 16.06.2014
 *      Author: Marc Nemes
 */

#ifndef EFFECTIVEDIAMETER_H_
#define EFFECTIVEDIAMETER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class EffectiveDiameter {

public:
	/**
	 * determines the effective diameter of a given graph
	 * the effective diameter is the smallest distance so that 90% percent of the nodes are connected to all other nodes within that distance
	 * @param G the given graph
	 * @return the effective diameter of the graph
	 */
	static count effectiveDiameter(const Graph& G);

};

} /* namespace NetworKit */

#endif /* ECCENTRICITY_H_ */
