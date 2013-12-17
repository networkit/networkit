/*
 * ClusteringCoefficient.h
 *
 *  Created on: 08.04.2013
 *      Author: cls
 */

#ifndef CLUSTERINGCOEFFICIENT_H_
#define CLUSTERINGCOEFFICIENT_H_

#include "../graph/Graph.h"

namespace NetworKit {

// TODO: is class necessary?
class ClusteringCoefficient {

public:

	ClusteringCoefficient();

	virtual ~ClusteringCoefficient();

	/**
	 * This calculates the average local clustering coefficient
	 * $$c(G) := \frac{1}{n} \sum_{u \in V} c(u)$$
	 *
	 * where $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$
	 */
	virtual double calculate(Graph& G);
};

} /* namespace NetworKit */
#endif /* CLUSTERINGCOEFFICIENT_H_ */
