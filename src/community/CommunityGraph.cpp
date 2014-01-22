/*
 * CommunityGraph.cpp
 *
 *  Created on: 22.08.2013
 *      Author: cls
 */

#include "CommunityGraph.h"

namespace NetworKit {

CommunityGraph::CommunityGraph() {
}

CommunityGraph::~CommunityGraph() {
}

void CommunityGraph::run(const Graph& G, const Clustering& zeta) {

	Gcom = Graph(0);
	Gcom.markAsWeighted(); // Gcon will be a weighted graph


	communityToSuperNode.clear();
	communityToSuperNode.resize(zeta.upperBound(), none); // there is one supernode for each cluster

	DEBUG("map cluster -> supernode");
	// populate map cluster -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		if (communityToSuperNode[c] == none) {
			communityToSuperNode[c] = Gcom.addNode();
		}
	});

	DEBUG("Gcom number of nodes: " , Gcom.numberOfNodes());


	index z = G.upperNodeIdBound();
	std::vector<node> nodeToSuperNode(z);

	DEBUG("node -> supernode");
	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = communityToSuperNode[zeta.clusterOf(v)];
	});


	DEBUG("create edges");

	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.forWeightedEdges([&](node u, node v, edgeweight ew) {
		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		if (Gcom.hasEdge(su, sv)) {
			Gcom.setWeight(su, su, Gcom.weight(su, su) + ew);
		} else {
			Gcom.addEdge(su, sv, ew);
		}
	});

}


Graph CommunityGraph::getGraph() {
	return Gcom;
}

std::map<index, node> CommunityGraph::getCommunityToNodeMap() {
	std::map<index, node> communityToNodeMap;
	for (index com = 0; com < communityToSuperNode.size(); ++com) {
		if (communityToSuperNode[com] != none) {
			communityToNodeMap[com] = communityToSuperNode[com];
		}
	}
	return communityToNodeMap;
}


std::map<node, index> CommunityGraph::getNodeToCommunityMap() {
	std::map<node, index> nodeToCommunityMap;
	for (index com = 0; com < communityToSuperNode.size(); ++com) {
		if (communityToSuperNode[com] != none) {
			node s = communityToSuperNode[com];
			nodeToCommunityMap[s] = com;
		}
	}
	return nodeToCommunityMap;
}

} /* namespace NetworKit */
