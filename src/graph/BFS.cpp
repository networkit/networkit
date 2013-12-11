/*
 * BFS.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include "BFS.h"

namespace NetworKit {


std::vector<count> BFS::run(const Graph& g, node source) const {
	count infDist = std::numeric_limits<count>::max();
	count n = g.numberOfNodes();
	std::vector<count> distances(n, infDist);
	std::queue<node> q;

	distances[source] = 0;
	q.push(source);

	while (! q.empty()) {
		node current = q.front();
		q.pop();
		TRACE("current node in BFS: " << current);

		// insert untouched neighbors into queue
		g.forNeighborsOf(current, [&](node neighbor) {
			if (distances[neighbor] == infDist) {
				q.push(neighbor);
				distances[neighbor] = distances[current] + 1;
			}
		});
	}

	return distances;
}

} /* namespace NetworKit */
