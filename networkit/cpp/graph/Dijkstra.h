/*
 * Dijkstra.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "Graph.h"
#include "SSSP.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dijkstra's SSSP algorithm.
 */
class Dijkstra : public SSSP {

friend class DynDijkstra;
friend class DynDijkstra2;

public:

	/**
	 * Creates the Dijkstra class for @a G and the source node @a source.
	 *
	 * @param G The graph.
	 * @param source The source node.
	 * @param storePaths	store paths and number of paths?
	 */
	Dijkstra(const Graph& G, node source, bool storePaths=true, bool storeStack=false, node target = none);

	/**
	 * Performs the Dijkstra SSSP algorithm on the graph given in the constructor.
	 */
	virtual void run();
};

} /* namespace NetworKit */
#endif /* DIJKSTRA_H_ */
