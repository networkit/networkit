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


// TODO: adapt BFS to standard interface: pass G and source via constructor, have a void run method, return results in appropriate getter methods.
class BFS {
public:
	BFS(const Graph& G, node source); 
	virtual ~BFS();

	/**
	 * Breadth-first search from @a source.
	 * @return Vector of unweighted distances from node @a source, i.e. the
	 * length (number of edges) of the shortest path from @a source to any other vertex.
	 */
	virtual void run();

	virtual std::vector<count> getDistances() const;

	virtual std::vector<node> getPath(node t, bool forward=true) const;

private: 
	const Graph& G;
	const node source;
	std::vector<count> distances;
	std::vector<node> previous;

};

} /* namespace NetworKit */
#endif /* BFS_H_ */
