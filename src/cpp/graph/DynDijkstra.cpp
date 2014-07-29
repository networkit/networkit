/*
 * DynDijkstra.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include "DynDijkstra.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include <queue>


namespace NetworKit {

DynDijkstra::DynDijkstra(const Graph& G, node s) : Dijkstra(G, s), color(G.upperNodeIdBound(), WHITE) {

}


void DynDijkstra::update(const std::vector<GraphEvent>& batch) {
	// priority queue with distance-node pairs
	Aux::PrioQueue<edgeweight, node> Q(G.upperNodeIdBound());
	// queue with all visited nodes
	std::queue<node> visited;
	// if u has a new shortest path going through v, it updates the distance of u
	// and inserts u in the priority queue (or updates its priority, if already in Q)
	auto updateQueue = [&](node u, node v, edgeweight w) {
		if (distances[u] >= distances[v]+w) {
			distances[u] = distances[v]+w;
			if (color[u] == WHITE) {
				Q.insert(distances[u], u);
				color[u] = BLACK;
			}	else {
				Q.decreaseKey(distances[u], u);
			}
		}
	};

	for (GraphEvent edge : batch) {
		if (edge.type!=GraphEvent::EDGE_ADDITION) //TODO: consider also weight decrease operations
			throw std::runtime_error("Graph update not allowed");
		updateQueue(edge.u, edge.v, edge.w);
		updateQueue(edge.v, edge.u, edge.w);
	}

	while(Q.size() != 0) {
		node current = Q.extractMin().second;
		visited.push(current);
		previous[current].clear();
		npaths[current] = 0;
		G.forWeightedNeighborsOf(current, [&](node z, edgeweight w){
			//z is a predecessor of current node
			if (distances[current] == distances[z]+w) {
				previous[current].push_back(z);
				npaths[current] += npaths[z];
			}
			//check whether curent node is a predecessor of z
			else {
				updateQueue(z, current, w);
			}
		});
	}

	// reset colors
	while(!visited.empty()) {
		node w = visited.front();
		visited.pop();
		color[w] = WHITE;
	}

}

} /* namespace NetworKit */
