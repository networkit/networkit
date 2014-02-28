/*
 * CoreDecomposition.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß, Christian Staudt
 */

#ifndef COREDECOMPOSITION_H_
#define COREDECOMPOSITION_H_

#include <vector>
#include <fstream>
#include <string>
#include <list>
#include "../graph/Graph.h"
#include "../auxiliary/ShellList.h"

namespace NetworKit {

/**
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition {

public:

	CoreDecomposition(const Graph& G);

	/**
	* Perform k-core decomposition of graph passed in constructor.
	*/
	void run();

	/**
	* @return vector or core numbers, indexed by node.
	*/
	std::vector<index> coreNumbers() const;

	/** @return core number of node @a v */
	index coreNumber(node v) const;

	/** @return the k-cores as sets of nodes, indexed by k */
	std::vector<std::set<node> > cores() const;

	/** @return the k-shells as sets of nodes, indexed by k */
	std::vector<std::set<node> > shells() const;


private:

	const Graph& G;
	std::vector<index> coreness;
	index maxCoreNumber; // maximum core number of any node in the graph
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_H_ */
