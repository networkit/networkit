/*
 * BFS.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef BFS_H_
#define BFS_H_

#include "Graph.h"

namespace NetworKit {

class BFS {
public:
	BFS() = default;
	virtual ~BFS() = default;

	/**
	 * Breadth-first search from @a source.
	 * @return Vector of unweighted distances from node @a source, i.e. the
	 * length (number of edges) of the shortest path from @a source to any other vertex.
	 */
	virtual std::vector<count> run(const Graph& g, node source) const;
};

} /* namespace NetworKit */
#endif /* BFS_H_ */
