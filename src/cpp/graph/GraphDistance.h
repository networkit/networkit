/*
 * GraphDistance.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef GRAPHDISTANCE_H_
#define GRAPHDISTANCE_H_

#include "Graph.h"
#include "Dijkstra.h"
#include "BFS.h"

namespace NetworKit {

// TODO: inherit from NodeDistance
class GraphDistance {
public:
	GraphDistance();
	virtual ~GraphDistance();

	/**
	 * @return Distance between @a u and @a v, i.e., the length of the shortest path
	 * between the two. Zero if u = v, maximal possible value if no path exists.
	 */
	virtual edgeweight weightedDistance(const Graph& g, node u, node v) const;

	/**
	 * @return Number of edges on shortest unweighted path between @a u and @a v.
	 */
	virtual count unweightedDistance(const Graph& g, node u, node v) const;
};

} /* namespace NetworKit */
#endif /* GRAPHDISTANCE_H_ */
