/*
 * DynBFS.cpp
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#include "DynBFS.h"
#include <queue>

namespace NetworKit {

DynBFS::DynBFS(const Graph& G, node s) : DynSSSP(G, s), BFS(G, s), color(G.upperNodeIdBound(), WHITE) {

}

void DynBFS::update(const std::vector<GraphEvent>& batch) {
	count maxLevels; // TODO: find out
	std::vector<std::queue<node> > queues(maxLevels);

	// insert nodes from the batch whose distance has changed (affected nodes) into the queues
	for (GraphEvent edge : batch) {
		if (edge.type!=GraphEvent::EDGE_ADDITION || edge.w!=1.0)
			throw std::runtime_error("wrong update");
		if (distances[edge.u] >= distances[edge.v]+1) {
			queues[distances[edge.v]+1].push(u);
		} else if (distances[edge.v] >= distances[edge.u]+1) {
			queues[distances[edge.u]+1].push(v);
		}
	}

	// extract nodes from the queues and scan incident edges
	std::queue<node> visited(G.upperNodeIdBound());
	int m = 1;
	while (m < maxLevels) {
		while (!queues[m].empty()) {
			node w = queues[m].pop();
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
					queues[m].push(z);
				}
			}
		}
		m = m+1;
	}

	// reset colors
	while(!visited.empty()) {
		node w = visited.pop();
		color[w] = WHITE;
	}
}

} /* namespace NetworKit */
