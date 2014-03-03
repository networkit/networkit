/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include "Dijkstra.h"

namespace NetworKit {

Dijkstra::Dijkstra() {

}

Dijkstra::~Dijkstra() {

}

std::vector<edgeweight> Dijkstra::run(const Graph& g, node source) {
	auto relax([&](node u, node v, edgeweight w, std::vector<edgeweight>& distances,
			Aux::PrioQueue<edgeweight, node>& pq)
	{
		if (distances[v] > distances[u] + g.weight(u, v)) {
			distances[v] = distances[u] + g.weight(u, v);
			pq.decreaseKey(distances[v], v);
		}
	});

	// init distances
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count n = g.numberOfNodes();
	std::vector<edgeweight> distances(n, infDist);
	distances[source] = 0;
	std::set<index> unsettled;

	// priority queue with distance-node pairs
	Aux::PrioQueue<edgeweight, node> pq(distances);

	while (pq.size() > 0) {
		DEBUG("pq size: ", pq.size());

		node current = pq.extractMin().second;
		DEBUG("pq size: ", pq.size());
		TRACE("current node in Dijkstra: " , current);

		g.forWeightedEdgesOf(current, [&](node current, node v, edgeweight w) {
			relax(current, v, w, distances, pq);
		});
	}

	return distances;
}

} /* namespace NetworKit */
