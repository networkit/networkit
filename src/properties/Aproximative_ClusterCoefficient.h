/*
 * Aproximative_ClusterCoefficient.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#ifndef APROXIMATIVE_CLUSTERCOEFFICIENT_H_
#define APROXIMATIVE_CLUSTERCOEFFICIENT_H_

#include <vector>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Computes k-core decomposition of a graph.
 */
class Aproximative_ClusterCoefficient {
public:
	Aproximative_ClusterCoefficient();
	virtual ~Aproximative_ClusterCoefficient();

	/**
	 * @return k-core decomposition of graph @a G.
	 */
	float run(const Graph& G,int k);
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */
