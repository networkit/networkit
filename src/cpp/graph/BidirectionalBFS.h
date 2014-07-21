/*
 * BidirectionalBFS.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Maximilian Vogel
 */

#ifndef BIDIRECTIONALBFS_H_
#define BIDIRECTIOANLBFS_H_

#include "Graph.h"
#include "BFS.h"

namespace NetworKit {

/**
 * @ingroup graph
 * The BFS class is used to do a breadth-first search on a Graph from a given source node.
 */
class BidirectionalBFS {
public:
	/**
	 * Constructs the BFS class for @a G and source node @a source.
	 *
	 * @param G The graph.
	 * @param source The source node of the breadth-first search.
	 */
	BidirectionalBFS(const Graph&); 

	/**
	 * Breadth-first search from @a source.
	 * @return Vector of unweighted distances from node @a source, i.e. the
	 * length (number of edges) of the shortest path from @a source to any other node.
	 */
	virtual void run(node s, node t);

private:
	const Graph& G;
	BFS forward, backward;

};

} /* namespace NetworKit */
#endif /* BFS_H_ */
