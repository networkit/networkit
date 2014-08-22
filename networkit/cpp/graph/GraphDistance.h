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
/**
 * @ingroup graph
 */
class GraphDistance {
public:

	/** Default destructor */
	virtual ~GraphDistance() = default;

	/**
	 * Returns the distance between @a u and @a v in Graph @a g i.e., the length of the shortest path
	 * between the two. Zero if u = v, maximal possible value if no path exists.
	 *
	 * @param g The graph.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @return The distance between @a u and @a v.
	 */
	virtual edgeweight weightedDistance(const Graph& g, node u, node v) const;

	/**
	 * Returns the number of edges on shortest unweighted path between @a u and @a v in Graph @a g.
	 *
	 * @param g The graph.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @return The number of edges between @a u and @a v.
	 */
	virtual count unweightedDistance(const Graph& g, node u, node v) const;
};

} /* namespace NetworKit */
#endif /* GRAPHDISTANCE_H_ */
