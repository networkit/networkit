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
 * Dijkstra's SSSP algorithm.
 */
class Dijkstra : public SSSP {

public:

	Dijkstra(const Graph& G, node source);
	
	virtual void run();



};

} /* namespace NetworKit */
#endif /* DIJKSTRA_H_ */
