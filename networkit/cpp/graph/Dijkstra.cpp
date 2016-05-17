/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#include "Dijkstra.h"

#include <algorithm>

namespace NetworKit {

Dijkstra::Dijkstra(const Graph& G, node source, bool storePaths, bool storeStack, node target) : SSSP(G, source, storePaths, storeStack, target) {

}

void Dijkstra::run() {

	TRACE("initializing Dijkstra data structures");
	// init distances
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	distances.clear();
	distances.resize(G.upperNodeIdBound(), infDist);
	if (storePaths) {
		previous.clear();
		previous.resize(G.upperNodeIdBound());
		npaths.clear();
		npaths.resize(G.upperNodeIdBound(), 0);
		npaths[source] = 1;
	}

	if (storeStack) {
		std::vector<node> empty;
		std::swap(stack, empty);
	}
	distances[source] = 0;
	// priority queue with distance-node pairs
	distances[source] = 0;
	Aux::PrioQueue<edgeweight, node> pq(distances);


	auto relax([&](node u, node v, edgeweight w) {
		if (distances[v] > distances[u] + w) {
			distances[v] = distances[u] + w;
			if (storePaths) {
				previous[v] = {u}; // new predecessor on shortest path
				npaths[v] = npaths[u];
			}
			TRACE("Decreasing key of ", v);
			TRACE("pq size: ", pq.size());
			pq.decreaseKey(distances[v], v);
			TRACE("pq size: ", pq.size());
		} else if (storePaths && (distances[v] == distances[u] + w)) {
			previous[v].push_back(u); 	// additional predecessor
			npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
		}
	});

	bool breakWhenFound = (target != none);
	TRACE("traversing graph");
	while (pq.size() > 0) {
		TRACE("pq size: ", pq.size());
		node current = pq.extractMin().second;
		TRACE("current node in Dijkstra: " , current);
		TRACE("pq size: ", pq.size());
		if (breakWhenFound && target == current) {
			break;
		}

		if (storeStack) {
			stack.push_back(current);
		}

		G.forEdgesOf(current, relax);
	}

}


} /* namespace NetworKit */
