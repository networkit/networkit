/*
 * Subgraph.cpp
 *
 *  Created on: Aug 26, 2014
 *      Author: Christian Staudt
 */

#include "Subgraph.h"

namespace NetworKit {

Graph Subgraph::fromNodes(const Graph& G, const std::unordered_set<node>& nodes) {

	Graph S(G.upperNodeIdBound(), G.isWeighted(), G.isDirected());
	// delete all nodes that are not in the node set
	S.forNodes([&](node u) {
		if (nodes.find(u) == nodes.end()) {
			S.removeNode(u);
		}
	});

	G.forEdges([&](node u, node v) {
		// if both end nodes are in the node set
		if (nodes.find(u) != nodes.end() && nodes.find(v) != nodes.end()) {
			S.addEdge(u, v, G.weight(u, v));
		}
	});

	return S;
}

} /* namespace NetworKit */
