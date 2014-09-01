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

namespace NetworKit {

/**
 * @ingroup properties
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition {

public:

	/**
	 * Create CoreDecomposition class for graph @a G.
	 *
	 * @param G The graph.
	 */
	CoreDecomposition(const Graph& G);

	/**
	* Perform k-core decomposition of graph passed in constructor.
	*/
	void run();

	/**
	 * Get vector of core numbers, indexed by node.
	 *
	 * @return Vector of core numbers, indexed by node.
	 */
	std::vector<index> coreNumbers() const;

	/**
	 * Get core number of node @a v.
	 *
	 * @param v A node.
	 * @return Core number of node @a v.
	 */
	index coreNumber(node v) const;

	/**
	 * Get the k-cores as sets of nodes, indexed by k.
	 *
	 * @return the k-cores as sets of nodes, indexed by k.
	 */
	std::vector<std::set<node> > cores() const;

	/**
	 * Get the k-shells as sets of nodes, indexed by k.
	 *
	 * @return the k-shells as sets of nodes, indexed by k
	 */
	std::vector<std::set<node> > shells() const;

	/**
	 * Get maximum core number.
	 *
	 * @return The maximum core number
	 */
	index maxCoreNumber() const;


private:

	const Graph& G;
	std::vector<index> coreness;
	index maxCore; // maximum core number of any node in the graph
	bool ran; // whether algorithm has been run
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */
