/*
 * BFS.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include "BFS.h"

namespace NetworKit {

BFS::BFS(const Graph& G, node source) : SSSP(G, source) {
}


void BFS::run() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	distances.clear();
	distances.resize(z, infDist);
	previous.clear();
	previous.resize(z);
	npaths.clear();
	npaths.resize(z, 0);

	std::queue<node> q;

	distances[source] = 0;
	npaths[source] = 1;
	q.push(source);

	while (! q.empty()) {
		node u = q.front();
		q.pop();
		// TRACE("current node in BFS: " , u);
//		TRACE(distances);

		// insert untouched neighbors into queue
		G.forNeighborsOf(u, [&](node v) {
			// TRACE("scanning neighbor ", v);

			if (distances[v] == infDist) {
				q.push(v);
				previous[v] = {u};
				distances[v] = distances[u] + 1;
				npaths[v] = npaths[u];
			} else if (distances[v] == distances[u] + 1) {
				previous[v].push_back(u); 	// additional predecessor
				npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
			}
		});
	}
}

} /* namespace NetworKit */
