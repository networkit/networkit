/*
 * Dijkstra.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "Graph.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * Dijkstra's SSSP algorithm.
 */
class Dijkstra {
protected:


public:
	Dijkstra();
	virtual ~Dijkstra();

	/**
	 * Dijkstra's SSSP algorithm.
	 * @return Vector of weighted distances from node @a source, i.e. the
	 * length of the shortest path from @a source to any other vertex.
	 */
	virtual std::vector<edgeweight> run(const Graph& g, node source);
};

} /* namespace NetworKit */
#endif /* DIJKSTRA_H_ */
