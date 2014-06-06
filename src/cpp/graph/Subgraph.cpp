/*
 * Subgraph.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: forigem
 */

#include "Subgraph.h"

namespace NetworKit {

Subgraph::Subgraph() {
	// TODO Auto-generated constructor stub

}

Subgraph::~Subgraph() {
	// TODO Auto-generated destructor stub
}

Graph Subgraph::fromNodes(const Graph& G, const std::unordered_set<node>& nodes) {

	Graph returnGraph;
	std::unordered_map<node,node> setNodesToSubgraphNodes;

	int i = 0;
	for (node setNode: nodes) {
		i++;
		node addedNode = returnGraph.addNode();
		setNodesToSubgraphNodes[setNode] = addedNode;
	}


	G.forEdges([&](node u, node v) {

		bool gotU = true;
		bool gotV = true;

		auto got = nodes.find(u);
		if (got == nodes.end())
			 gotU = false;

		if (gotU == true) {
			auto got = nodes.find(v);
			if (got == nodes.end())
				gotV = false;
		}

		if (gotV == true && gotU == true) {
			node uSubgraph = setNodesToSubgraphNodes[u];
			node vSubgraph = setNodesToSubgraphNodes[v];
			returnGraph.addEdge(uSubgraph, vSubgraph);
		}


	});
	return returnGraph;
}

} /* namespace NetworKit */
