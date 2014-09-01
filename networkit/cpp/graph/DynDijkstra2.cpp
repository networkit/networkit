/*
 * DynDijkstra2.cpp
 *
 *  Created on: 24.07.2014
 *      Author: ebergamini
 */

#include "DynDijkstra2.h"
#include "Dijkstra.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include <queue>


namespace NetworKit {

DynDijkstra2::DynDijkstra2(const Graph& G, node source, bool storePredecessors) : DynSSSP(G, source, storePredecessors),
color(G.upperNodeIdBound(), WHITE),
N (G.upperNodeIdBound(), Aux::PrioQueue<double, node>(G.upperNodeIdBound())){
}


void DynDijkstra2::run(node t) {
	INFO("Running dijkstra to initialize the structures");
	if (t != none) {
		throw std::runtime_error("Invalid argument: DynDijkstra2 doesn't work with a target node.");
	}
	Dijkstra dij(G, source, true);
	dij.run();
	INFO("Finished running dijkstra");
	distances = dij.distances;
	npaths = dij.npaths;
	if (storePreds) {
		previous = dij.previous;
	}
	// insert all neighbors of each vertex u in its priority queue
	G.forNodes([&](node u){
		G.forNeighborsOf(u, [&](node v, edgeweight w){
			N[u].insert(distances[v]-w, v);
		});
	});
}

void DynDijkstra2::update(const std::vector<GraphEvent>& batch) {
	// priority queue with distance-node pairs
	Aux::PrioQueue<edgeweight, node> Q(G.upperNodeIdBound());
	// queue with all visited edges
	std::queue<std::pair<node, node>> visited;
	std::vector<count> old_paths = npaths;
	// if u has a new shortest path going through v, it updates the distance of u
	// and inserts u in the priority queue (or updates its priority, if already in Q)
	auto updateQueueAndPaths = [&](node u, node v, edgeweight w) {
		if (distances[u] > distances[v]+w) {
			DEBUG("Found shorter path going through ", v);
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
			DEBUG("Found new shortest paths (not shorter) going through ", v);
			// check whether v was already a predecessor of u (in that case, subtract the old contribution)
			bool old_pred = false;
			for (node p : previous[u]) {
				if (p == v)
					old_pred = true;
			}
			if (old_pred) {
				npaths[u] -= old_paths[v];
			} else {
				previous[u].push_back(v);
			}
			npaths[u] += npaths[v];
			if (color[u] == WHITE) {
				Q.insert(distances[u], u);
				color[u] = BLACK;
			}
		}
	};

	INFO("Inserting edges of the batch in the priority queues N");
	for (GraphEvent edge : batch) {
		if (edge.type!=GraphEvent::EDGE_ADDITION) //TODO: consider also weight decrease operations
			throw std::runtime_error("Graph update not allowed");
		// insert edge.u in edge.v's priority queue and vice versa (if you allow also weight decrease, this will need to be changed)
		N[edge.v].insert(-1* (distances[edge.u]-edge.w), edge.u);
		N[edge.u].insert(-1 * (distances[edge.v]-edge.w), edge.v);
		updateQueueAndPaths(edge.u, edge.v, edge.w);
		updateQueueAndPaths(edge.v, edge.u, edge.w);
	}
	INFO("Extracting vertices from the main priority queue");
	while(Q.size() != 0) {
		node current = Q.extractMin().second;
		DEBUG("Extracted node ", current);
		visited.push(std::pair<node, node>(current, current));
		std::pair<double, node> neighbor = N[current].extractMin();
		edgeweight priority = -1 * neighbor.first;
		while((N[current].size() > 0) && (priority >= distances[current])){
			DEBUG("Neighbor: ", neighbor.second);
			visited.push(std::pair<node, node>(current, neighbor.second));
			updateQueueAndPaths(neighbor.second, current, G.weight(current, neighbor.second));
			if (N[current].size() > 0) {
				neighbor = N[current].extractMin();
				priority = -1 * neighbor.first;
			}
		}

	}
	INFO("Resetting colors");
	// reset colors
	while(!visited.empty()) {
		std::pair<node, node> e = visited.front();
		visited.pop();
		color[e.first] = WHITE;
		if (e.first != e.second) {
			color[e.second] = WHITE;
			N[e.first].insert(-1*(distances[e.second]-G.weight(e.first, e.second)), e.second);
			//TODO shouldn't we also update distances in the priority queue of e.second? Because e.first has changed its distance...
		}
	}

}

} /* namespace NetworKit */
