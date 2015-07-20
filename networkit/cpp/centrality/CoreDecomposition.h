/*
 * CoreDecomposition.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiss, Christian Staudt
 */

#ifndef COREDECOMPOSITION_H_
#define COREDECOMPOSITION_H_

#include <vector>
#include <fstream>
#include <string>
#include <list>
#include "../graph/Graph.h"
#include "../centrality/Centrality.h"
#include "../structures/Partition.h"
#include "../structures/Cover.h"


namespace NetworKit {

/**
 * @ingroup centrality
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition : public NetworKit::Centrality  {

public:

	/**
	 * Create CoreDecomposition class for graph @a G. The graph may not contain self-loops.
	 *
	 * @param G The graph.
	 */
	CoreDecomposition(const Graph& G);

	/**
	* Perform k-core decomposition of graph passed in constructor.
	*/
	void run();

	/**
	 * Get the k-cores as a graph cover object.
	 *
	 * @return the k-cores as a Cover
	 */
	Cover cores() const;

	/**
	 * Get the k-shells as a partition object
	 *
	 * @return the k-shells as a Partition
	 */
	Partition shells() const;

	/**
	 * Get maximum core number.
	 *
	 * @return The maximum core number
	 */
	index maxCoreNumber() const;

private:

	index maxCore; // maximum core number of any node in the graph
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */
