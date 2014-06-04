/*
 * MSTFinder.h
 *
 *  Created on: 13.05.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MSTFINDER_H_
#define MSTFINDER_H_

#include <vector>
#include <utility>

#include "../graph/Graph.h"

namespace NetworKit {

typedef std::pair<node, node> Edge;

/**
 * Abstract base class for MST algorithms.
 */
class MSTFinder {
public:
	/** destructor */
	virtual ~MSTFinder() {}

	/**
	 * Computes an MST for each connected component of @a G.
	 * @return A vector of MSTs, each being a vector of edges (node pairs).
	 */
	virtual std::vector<std::vector<Edge>> run(const Graph &G) = 0;
};

} /* namespace NetworKit */

#endif /* MSTFINDER_H_ */
