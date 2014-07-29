/*
 * DynDijkstra2.cpp
 *
 *  Created on: 24.07.2014
 *      Author: ebergamini
 */

#include "DynDijkstra2.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include <queue>


namespace NetworKit {

DynDijkstra2::DynDijkstra2(const Graph& G, node s) : Dijkstra(G, s),
color(G.upperNodeIdBound(), WHITE),
N (G.upperNodeIdBound(), Aux::PrioQueue<edgeweight, node>(G.upperNodeIdBound())){

	// insert all neighbors of each vertex u in its priority queue
	G.forNodes([&](node u){
		G.forWeightedNeighborsOf(u, [&](node v, edgeweight w){
			N[u].insert(distances[v]-w, v);
		});
	});
}


void DynDijkstra2::update(const std::vector<GraphEvent>& batch) {
	// priority queue with distance-node pairs
/*	Aux::PrioQueue<edgeweight, node> Q(G.upperNodeIdBound());
	// queue with all visited edges
	std::queue<edge> visited;
	// if u has a new shortest path going through v, it updates the distance of u
	// and inserts u in the priority queue (or updates its priority, if already in Q)
	auto updateQueueAndPaths = [&](node u, node v, edgeweight w) {
		if (distances[u] > distances[v]+w) {
			distances[u] = distances[v]+w;
			previous[u] = {v};
			npaths[u] = npaths[v];
			if (color[u] == WHITE) {
				Q.insert(distances[u], u);
				color[u] = BLACK;
			}	else {
				Q.decreaseKey(distances[u], u);
			}
		} else if (distances[u] == distances[v]+w) {
			previous[u].push_back(v);
			npaths[u] += napths[v];
			if (color[u] == WHITE) {
				Q.insert(distances[u], u);
				color[u] = BLACK;
			}
		}
	};

	for (GraphEvent edge : batch) {
		if (edge.type!=GraphEvent::EDGE_ADDITION) //TODO: consider also weight decrease operations
			throw std::runtime_error("Graph update not allowed");
		// insert edge.u in edge.v's priority queue and vice versa
		N[edge.v].insert(distances[edge.u]-edge.w, edge.u);
		N[edge.u].insert(distances[edge.v]-edge.w, edge.v);
		updateQueue(edge.u, edge.v, edge.w);
		updateQueue(edge.v, edge.u, edge.w);
	}

	while(Q.size() != 0) {
		// NO!! here we need a max-based priority queue!!
		node current = Q.extractMin().second;
		neighbor = N[current].extractMin();
		while(neighbor.first >= distances[current]){

		}

	}

	// reset colors
	while(!visited.empty()) {
		edge e = visited.front();
		visited.pop();
		color[e.u] = WHITE;
		color[e.v] = WHITE;
	}*/

}

} /* namespace NetworKit */
