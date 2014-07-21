/*
 * BFS.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include <queue>

#include "BFS.h"
#include "../auxiliary/Log.h"

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

//	std::queue<node> q;

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

node BFS::settleNext() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	node current = q.front();
	q.pop();
	G.forNeighborsOf(current, [&](node v) {
		// TRACE("scanning neighbor ", v);
		if (distances[v] == infDist) {
			q.push(v);
			previous[v] = {current};
			distances[v] = distances[current] + 1;
			npaths[v] = npaths[current];
		} else if (distances[v] == distances[current] + 1) {
			previous[v].push_back(current); 	// additional predecessor
			npaths[v] += npaths[current]; 	// all the shortest paths to u are also shortest paths to v now
		}
	});
	return current;
}

void BFS::init(node s) {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	distances.clear();
	distances.resize(z, infDist);
	previous.clear();
	previous.resize(z);
	npaths.clear();
	npaths.resize(z, 0);

//	std::queue<node> q;

	distances[s] = 0;
	npaths[s] = 1;
	q.push(s);
}

node BFS::isFinished() {
	return q.empty();
}

node BFS::getCurrentMin() {
	return distances[q.front()];
}

bool BFS::wasVisited(node u) {
	return distances[u] < std::numeric_limits<edgeweight>::max();
}

edgeweight BFS::extractDistance(node u) {
	return distances[u];
}


void BFS::runUntil(node t) {
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
	bool found = false;
	while (! q.empty() && !found) {
		node u = q.front();
		if (t == u) {
			found = true;
			continue;
		}
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

void BFS::runBidirectional(node t) {
	TRACE("does this method even gets called?");
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	distances.clear();
	distances.resize(z, infDist);
	previous.clear();
	previous.resize(z);
	npaths.clear();
	npaths.resize(z, 0);

	std::vector<double> distancesBack(z, infDist);
	std::vector<std::vector<node>> previousBack(z);
	std::vector<count> npathsBack(z,0);

	std::queue<node> q;
	std::queue<node> qBack;

	distances[source] = 0;
	distancesBack[t] = 0;
	npathsBack[t] = 1;
	qBack.push(t);
	npaths[source] = 1;
	q.push(source);
	bool found = false;
	bool forwardNext = true;

	std::vector<node> junction;
	double minDist = infDist;
//	auto continueSearch = [&](double minDistance, double forwardMin, double backwardMin) {
//		if (distances[forwardMin] == infDist || distancesBack[backwardMin] == infDist)
//			return true;
//		return (minDistance > (distances[forwardMin] + distancesBack[backwardMin]));
//	};
	TRACE("starting bidirectional BFS");
	while (! q.empty() && !qBack.empty()) { //&& continueSearch(minDist,q.front(),qBack.front())) {
		bool path_found = false;
		node u;
		if (forwardNext) {
			u = q.front();
			TRACE("forward search with node ",u);
			if (t == u) {
				found = true;
				continue;
			}
			q.pop();
			// insert untouched neighbors into queue
			G.forNeighborsOf(u, [&](node v) {
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
			if (distancesBack[u] != infDist) path_found = true;
		} else {
			u = qBack.front();
			TRACE("backward search with node ",u);
			if (t == u) {
				found = true;
				continue;
			}
			qBack.pop();
			G.forNeighborsOf(u, [&](node v) {
				if (distancesBack[v] == infDist) {
					qBack.push(v);
					previousBack[v] = {u};
					distancesBack[v] = distances[u] + 1;
					npathsBack[v] = npathsBack[u];
				} else if (distancesBack[v] == distancesBack[u] + 1) {
					previousBack[v].push_back(u); 	// additional predecessor
					npathsBack[v] += npathsBack[u]; 	// all the shortest paths to u are also shortest paths to v now
				}
			});
			if (distances[u] != infDist) path_found = true;
		}
		if (path_found) {
			double distance = distances[u] + distancesBack[u];
			if (distance < minDist) {
				TRACE("new minimal path found");
				minDist = distance;
				junction = {u};
			} else if (distance == minDist) {
				TRACE("new node with same distance found");
				junction.push_back(u);
			}
		}
		forwardNext = !forwardNext;
	}


}



} /* namespace NetworKit */
