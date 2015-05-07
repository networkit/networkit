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

namespace NetworKit {

/**
 * @ingroup properties
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition : public NetworKit::Centrality  {

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

	index maxCore; // maximum core number of any node in the graph
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */
