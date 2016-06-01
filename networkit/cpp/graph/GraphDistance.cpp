/*
 * GraphDistance.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include "GraphDistance.h"

namespace NetworKit {

edgeweight GraphDistance::weightedDistance(const Graph& g, node u, node v) const {
	Dijkstra dijkstra(g, u); 
	dijkstra.run();
	std::vector<edgeweight> distances = dijkstra.getDistances();
	DEBUG("Called Dijkstra, distance between " , u , " and " , v , ": " , distances[v]);
	return distances[v];
}

count GraphDistance::unweightedDistance(const Graph& g, node u, node v) const {
	BFS bfs(g, u);
	DEBUG("running BFS");
	bfs.run();
	auto distances = bfs.getDistances();
	DEBUG("Called BFS, distance between " , u , " and " , v , ": " , distances[v]);
	return (count) distances[v];
}

} /* namespace NetworKit */
