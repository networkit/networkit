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

Graph Subgraph::fromNodes(Graph G,
		std::unordered_set<node> nodeMap) {

	Graph returnGraph;
	std::unordered_map<node,node> setNodesToSubgraphNodes;

	int i = 0;
	for (node setNode: nodeMap) {
		DEBUG(i);
		i++;
		node addedNode = returnGraph.addNode();
		setNodesToSubgraphNodes[setNode] = addedNode;
	}


	G.parallelForEdges([&](node u, node v) {

		bool gotU = true;
		bool gotV = true;

		auto got = nodeMap.find(u);
		if (got == nodeMap.end())
			 gotU = false;

		if (gotU == true) {
			auto got = nodeMap.find(v);
			if (got == nodeMap.end())
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
