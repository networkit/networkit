/*
 * GraphDistance.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include "GraphDistance.h"

namespace NetworKit {

GraphDistance::GraphDistance() {
	// TODO Auto-generated constructor stub

}

GraphDistance::~GraphDistance() {
	// TODO Auto-generated destructor stub
}

edgeweight GraphDistance::weightedDistance(const Graph& g, node u, node v) const {
	Dijkstra dijkstra;
	std::vector<edgeweight> distances = dijkstra.run(g, u);
	DEBUG("Called Dijkstra, distance between " << u << " and " << v << ": " << distances[v]);
	return distances[v];
}

count GraphDistance::unweightedDistance(const Graph& g, node u, node v) const {
	BFS bfs;
	std::vector<count> distances = bfs.run(g, u);
	DEBUG("Called BFS, distance between " << u << " and " << v << ": " << distances[v]);
	return distances[v];
}

} /* namespace NetworKit */
