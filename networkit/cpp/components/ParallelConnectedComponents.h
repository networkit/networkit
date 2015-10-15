/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#ifndef PARALLELCONNECTEDCOMPONENTS_H_
#define PARALLELCONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../graph/BFS.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup components
 * Determines the connected components of an undirected graph.
 */
class ParallelConnectedComponents : public Algorithm {
public:

	ParallelConnectedComponents(const Graph& G, bool coarsening = true);

	/**
	 * This method determines the connected components for the graph g.
	 */
	void runSequential();

	/**
	 * This method determines the connected components for the graph g.
	 */
	void run();

	/**
	 * This method returns the number of connected components.
	 */
	count numberOfComponents();

	/**
	 * This method returns the the component in which node query is situated.
	 *
	 * @param[in]	query	the node whose component is asked for
	 */
	count componentOfNode(node u);


	/**
	 * Return a Partition that represents the components
	 */
	Partition getPartition();


private:
	const Graph& G;
	Partition component;
	bool coarsening;
};

}


#endif /* PARALLELCONNECTEDCOMPONENTS_H_ */
