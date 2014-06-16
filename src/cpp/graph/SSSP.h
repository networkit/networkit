/*
 * SSSP.h
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#ifndef SSSP_H_
#define SSSP_H_

#include <set>
 
#include "Graph.h"

namespace NetworKit {

/** 
 * Abstract base class for single-source shortest path algorithms.
 */
class SSSP {

public:

	SSSP(const Graph& G, node s);

	virtual ~SSSP() = default;	

	virtual void run() = 0;

	/**
	 * return Vector of weighted distances from node @a source, i.e. the
 	 * length of the shortest path from @a source to any other vertex.
	 */
	virtual std::vector<edgeweight> getDistances() const;

	/**
	 * @param  t target node
	 * @return   distance from s to target node t
	 * 	 */
	virtual edgeweight distance(node t) const;

	/**
	 * @param  t target node
	 * @return   number of shortest paths between s and t
	 * 	 */
	virtual count numberOfPaths(node t) const;

	/**
	 * @param t target node
	 * @return predecessors of t on all shortest paths from source to t
	 */
	virtual std::vector<node> getPredecessors(node t) const;

	/**
	 * @return a shortest path from source node to target node @a t.
	 * Returns empty path if source and target are not connected.
	 */
	virtual std::vector<node> getPath(node t, bool forward=true) const;

	/**
	 * @return all shortest paths from source node to target node @a t.
	 * Returns empty set if source and target are not connected.
	 */
	virtual std::set<std::vector<node> > getPaths(node t, bool forward=true) const;

protected:

	const Graph& G;
	const node source;
	std::vector<edgeweight> distances;
	std::vector<std::vector<node> > previous; // predecessors on shortest path
	std::vector<count> npaths;
};

} /* namespace NetworKit */

#endif /* SSSP_H_ */
