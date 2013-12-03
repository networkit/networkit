/*
 * CommunityGraph.cpp
 *
 *  Created on: 22.08.2013
 *      Author: cls
 */

#include "CommunityGraph.h"

namespace NetworKit {

CommunityGraph::CommunityGraph() {
	// TODO Auto-generated constructor stub

}

CommunityGraph::~CommunityGraph() {
	// TODO Auto-generated destructor stub
}

std::pair<Graph, std::vector<node> > CommunityGraph::get(const Graph& G, const Clustering& zeta) {
	Graph Gcom(0); // empty graph
	Gcom.markAsWeighted(); // Gcon will be a weighted graph

	std::vector<node> communityToSuperNode(zeta.upperBound(), none); // there is one supernode for each cluster

	DEBUG("map cluster -> supernode");
	// populate map cluster -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		if (communityToSuperNode[c] == none) {
			communityToSuperNode[c] = Gcom.addNode();
		}
	});

	DEBUG("Gcom number of nodes: " << Gcom.numberOfNodes());


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

	return std::make_pair(Gcom, communityToSuperNode);

}

} /* namespace NetworKit */
