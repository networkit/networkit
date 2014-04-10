/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#include "Dijkstra.h"

namespace NetworKit {

Dijkstra::Dijkstra(const Graph& G, node source) : G(G), source(source) {

}

Dijkstra::~Dijkstra() {

}

void Dijkstra::run() {

	// init distances
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	distances.clear();
	distances.resize(G.upperNodeIdBound(), infDist);
	previous.clear();
	previous.resize(G.upperNodeIdBound(), none); 
	distances[source] = 0;

	// priority queue with distance-node pairs
	Aux::PrioQueue<edgeweight, node> pq(distances);


	auto relax([&](node u, node v, edgeweight w) {
		if (distances[v] > distances[u] + w) {
			distances[v] = distances[u] + w;
			previous[v] = u; // new predecessor on shortest path
			pq.decreaseKey(distances[v], v);
		}
	});

	while (pq.size() > 0) {
//		DEBUG("pq size: ", pq.size());

		node current = pq.extractMin().second;
//		DEBUG("pq size: ", pq.size());
//		TRACE("current node in Dijkstra: " , current);

		G.forWeightedEdgesOf(current, relax);
	}

}


std::vector<edgeweight> Dijkstra::getDistances() const {
	return distances;
}


std::vector<node> Dijkstra::getPath(node t, bool forward) const {
	std::vector<node> path;
	node v = t;
	while (v != source) {
		path.push_back(v);
		v = previous[t];
	}

	if (forward) {
		std::reverse(path.begin(), path.end());
	}
	return path;
}

} /* namespace NetworKit */
