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

Graph CommunityGraph::get(const Graph& G, const Clustering& zeta) {
	Graph Gcom(0); // empty graph
	Gcom.markAsWeighted(); // Gcon will be a weighted graph

	IndexMap<cluster, node> clusterToSuperNode(zeta.upperBound(), none); // there is one supernode for each cluster

	// populate map cluster -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		if (! clusterToSuperNode.hasBeenSet(c)) {
			clusterToSuperNode[c] = Gcom.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});


	count n = G.numberOfNodes();
	NodeMap<node> nodeToSuperNode(n);

	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = clusterToSuperNode[zeta.clusterOf(v)];
	});


	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.forWeightedEdges([&](node u, node v, edgeweight ew) {
		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			// add edge weight to supernode (self-loop) weight
			Gcom.setWeight(su, su, Gcom.weight(su, su) + ew);
		} else {
			// add edge weight to weight between two supernodes (or insert edge)
			Gcom.setWeight(su, sv, Gcom.weight(su, sv) + ew);
		}
	});

	return Gcom;

}

} /* namespace NetworKit */
