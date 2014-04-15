/*
 * BFS.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef BFS_H_
#define BFS_H_

#include "Graph.h"
#include "SSSP.h"

namespace NetworKit {


class BFS : public SSSP {
public:
	BFS(const Graph& G, node source); 

	/**
	 * Breadth-first search from @a source.
	 * @return Vector of unweighted distances from node @a source, i.e. the
	 * length (number of edges) of the shortest path from @a source to any other vertex.
	 */
	virtual void run();

};

} /* namespace NetworKit */
#endif /* BFS_H_ */
