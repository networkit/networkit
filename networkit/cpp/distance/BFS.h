/*
 * BFS.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef BFS_H_
#define BFS_H_

#include "../graph/Graph.h"
#include "SSSP.h"

namespace NetworKit {

/**
 * @ingroup distance
 * The BFS class is used to do a breadth-first search on a Graph from a given source node.
 */
class BFS : public SSSP {

friend class DynBFS;

public:
	/**
	 * Constructs the BFS class for @a G and source node @a source.
	 *
	 * @param G The graph
	 * @param source The source node of the breadth-first search
	 * @param storePaths Paths are reconstructable and the number of paths is stored.
	 * @param storeNodesSortedByDistance Store a vector of nodes ordered in increasing distance from the source.
	 * @param target The target node.
	 */
	BFS(const Graph& G, node source, bool storePaths=true, bool storeNodesSortedByDistance=false, node target = none);

	/**
	 * Breadth-first search from @a source.
	 * @return Vector of unweighted distances from node @a source, i.e. the
	 * length (number of edges) of the shortest path from @a source to any other node.
	 */
	virtual void run();

};

} /* namespace NetworKit */
#endif /* BFS_H_ */
