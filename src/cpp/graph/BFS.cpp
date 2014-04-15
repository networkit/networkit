/*
 * BFS.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include "BFS.h"

namespace NetworKit {

BFS::BFS(const Graph& G, node source) : G(G), source(source) {
}

BFS::~BFS() {

}

void BFS::run() {
	count z = G.upperNodeIdBound();
	distances.clear();
	distances.resize(z, none);
	previous.clear();
	previous.resize(z);
	std::queue<node> q;

	distances[source] = 0;
	q.push(source);
	previous[source] = source;

	while (! q.empty()) {
		node current = q.front();
		q.pop();
		TRACE("current node in BFS: " , current);

		// insert untouched neighbors into queue
		G.forNeighborsOf(current, [&](node neighbor) {
			if (distances[neighbor] == none) {
				q.push(neighbor);
				previous[neighbor] = current;
				distances[neighbor] = distances[current] + 1;
			}
		});
	}
}

std::vector<count> BFS::getDistances() const {
	return distances;
}

std::vector<node> BFS::getPath(node t, bool forward) const {
	std::vector<node> path;
	if (previous[t] == none) { // t is not reachable from source
		WARN("there is no path from ", source, " to ", t);
		return path;
	}
	node v = t;
	while (v != source) {
		path.push_back(v);
		v = previous[v];
	}
	path.push_back(v); // appends source node, probably not necessary.

	if (forward) {
		std::reverse(path.begin(), path.end());
	}
	return path;
}

} /* namespace NetworKit */
