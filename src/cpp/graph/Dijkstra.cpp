/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#include "Dijkstra.h"

#include <algorithm>

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

edgeweight Dijkstra::distance(node t) const {
	return distances[t];	
}

count Dijkstra::numberOfPaths(node t) const {
	return npaths[t];	
}


std::vector<node> Dijkstra::getPath(node t, bool forward) const {
	std::vector<node> path;
	if (previous[t].empty()) { // t is not reachable from source
		WARN("there is no path from ", source, " to ", t);
		return path;
	}
	node v = t;
	while (v != source) {
		path.push_back(v);
		v = previous[t].front();
	}

	if (forward) {
		std::reverse(path.begin(), path.end());
	}
	return path;
}


std::set<std::vector<node> > Dijkstra::getPaths(node t, bool forward) const {
	throw std::runtime_error("FIXME: correct implementation needed");

	std::set<std::vector<node> > paths;
	if (previous[t].empty()) { // t is not reachable from source
		WARN("there is no path from ", source, " to ", t);
		return paths;
	}


	std::function<std::set<std::vector<node> > (std::vector<node>& prefix, node v) > trace = [&](std::vector<node>& prefix, node v) {
		// base case

		prefix.push_back(v);
		std::set<std::vector<node> > paths;
		paths.insert(prefix);
		for (node u : previous[v]) {
			auto returned = trace(prefix, u);
			paths.insert(returned.begin(), returned.end());
		}
		return paths;
	};

	std::vector<node> emptyPath;
	auto thePaths = trace(emptyPath, t);

	if (forward) {
		for (auto path : thePaths) {
			std::reverse(path.begin(), path.end());
		}
	}
	return thePaths;
}

} /* namespace NetworKit */
