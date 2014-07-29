/*
 * DynBFS.cpp
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#include "DynBFS.h"
#include "../auxiliary/Log.h"
#include <queue>


namespace NetworKit {

DynBFS::DynBFS(const Graph& G, node s) : BFS(G, s), color(G.upperNodeIdBound(), WHITE) {
}

void DynBFS::init() {
	run();
	maxDistance = 0;
	G.forNodes([&] (node n){
		if (distances[n] > maxDistance)
			maxDistance = distances[n];
	});
	maxDistance++;
}

void DynBFS::update(const std::vector<GraphEvent>& batch) {
	std::vector<std::queue<node> > queues(maxDistance);

	// insert nodes from the batch whose distance has changed (affected nodes) into the queues
	for (GraphEvent edge : batch) {
		if (edge.type!=GraphEvent::EDGE_ADDITION || edge.w!=1.0)
			throw std::runtime_error("Graph update not allowed");
		if (distances[edge.u] >= distances[edge.v]+1) {
			queues[distances[edge.v]+1].push(edge.u);
		} else if (distances[edge.v] >= distances[edge.u]+1) {
			queues[distances[edge.u]+1].push(edge.v);
		}
	}

	// extract nodes from the queues and scan incident edges
	std::queue<node> visited;
	count m = 1;
	while (m < maxDistance) {
		DEBUG("m = ", m);
		while (!queues[m].empty()) {
			node w = queues[m].front();
			DEBUG("node ", w);
			queues[m].pop();
			if (color[w] == BLACK) {
				continue;
			}
			visited.push(w);
			color[w] = BLACK;
			distances[w] = m;
			previous[w].clear();
			npaths[w] = 0;
			G.forNeighborsOf(w, [&](node z) {
				//z is a predecessor for w
				if (distances[w] == distances[z]+1) {
					previous[w].push_back(z);
					npaths[w] += npaths[z];
				}
				//w is a predecessor for z
				else if (color[z] == WHITE && distances[z] >= distances[w]+1 ) {
					color[z] = GRAY;
					queues[m+1].push(z);
				}
			});
		}
		m = m+1;
	}

	// reset colors
	while(!visited.empty()) {
		node w = visited.front();
		visited.pop();
		color[w] = WHITE;
	}
}

} /* namespace NetworKit */
