/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#include "Dijkstra.h"

#include <algorithm>

namespace NetworKit {

Dijkstra::Dijkstra(const Graph& G, node source) : SSSP(G, source) {

}




void Dijkstra::run() {

	DEBUG("initializing Dijkstra data structures");
	// init distances
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	distances.clear();
	distances.resize(G.upperNodeIdBound(), infDist);
	previous.clear();
	previous.resize(G.upperNodeIdBound()); 
	distances[source] = 0;
	npaths.clear();
	npaths.resize(G.upperNodeIdBound(), 0);
	npaths[source] = 1;

	// priority queue with distance-node pairs
	Aux::PrioQueue<edgeweight, node> pq(distances);


	auto relax([&](node u, node v, edgeweight w) {
		if (distances[v] > distances[u] + w) {
			distances[v] = distances[u] + w;
			previous[v] = {u}; // new predecessor on shortest path
			npaths[v] = npaths[u];
			pq.decreaseKey(distances[v], v);
		} else if (distances[v] == distances[u] + w) {
			previous[v].push_back(u); 	// additional predecessor
			npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
		}
	});


	DEBUG("traversing graph");
	while (pq.size() > 0) {
//		DEBUG("pq size: ", pq.size());

		node current = pq.extractMin().second;
//		DEBUG("pq size: ", pq.size());
//		TRACE("current node in Dijkstra: " , current);

		G.forEdgesOf(current, relax);
	}

}

} /* namespace NetworKit */
