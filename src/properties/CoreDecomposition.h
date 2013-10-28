/*
 * CoreDecomposition.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#ifndef COREDECOMPOSITION_H_
#define COREDECOMPOSITION_H_

#include <vector>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition {
public:
	CoreDecomposition();
	virtual ~CoreDecomposition();

	/**
	 * @return k-core decomposition of graph @a G.
	 */
	std::vector<count> run(const Graph& G);
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */
