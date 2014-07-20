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
	 * approximates the effective diameter of a given graph
	 * the effective diameter is the smallest distance so that 90% percent of the nodes are connected to all other nodes within that distance
	 * @param G the given graph
	 * @param ratio the ratio of nodes that should to be connected (0,1]
	 * @param k the number of parallel approximations to get a more robust result
	 * @param r the amount of bits that are added to the length of the bitmask to improve the accuracy
	 * @param l the number of iterations that a node has to iterate without connecting to new neighbors before it is no longer considered
	 * @return the approximated effective diameter of the graph
	 */
	static count effectiveDiameter(const Graph& G, const double ratio, const count k, const count r);
	static count effectiveDiameter(const Graph& G, const double ratio);
	static count effectiveDiameter(const Graph& G);

	static count effectiveDiameterExact(const Graph& G);

};

} /* namespace NetworKit */

#endif /* ECCENTRICITY_H_ */
