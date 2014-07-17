/*
 * DynBFS.cpp
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#include "DynBFS.h"
#include <queue>

namespace NetworKit {

DynBFS::DynBFS(const Graph& G, node s) : DynSSP(G, s), BFS(G, s), color(G.upperNodeIdBound(), WHITE) {

}

void DynBFS::update(const std::vector<GraphEvent>& batch) {
	count maxLevels; // TODO: find out
	std::vector<std::queue<node> > queues(maxLevels);

	// insert nodes from the batch whose distance has changed (affected nodes) into the queues

	// extract nodes from the queues and scan incident edges


	// reset colors
}

} /* namespace NetworKit */
