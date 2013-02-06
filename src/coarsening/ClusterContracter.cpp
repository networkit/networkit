/*
 * ClusterContracter.cpp
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#include "ClusterContracter.h"

namespace EnsembleClustering {

ClusterContracter::ClusterContracter() {
	// TODO Auto-generated constructor stub

}

ClusterContracter::~ClusterContracter() {
	// TODO Auto-generated destructor stub
}

std::pair<Graph, NodeMap<node> > ClusterContracter::run(Graph& G, Clustering& zeta) {

	Graph Gcon(0); // empty graph

	IndexMap<cluster, node> clusterToSuperNode(zeta.upperBound(), 0); // there is one supernode for each cluster

	// populate map cluster -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		if (! clusterToSuperNode.hasBeenSet(c)) {
			node sv = Gcon.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
			clusterToSuperNode[c] = sv;
		}
	});


	int64_t n = G.numberOfNodes();
	NodeMap<node> nodeToSuperNode(n);

	// set entries node -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		nodeToSuperNode[v] = clusterToSuperNode[c];
	});


	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.forEdges([&](node u, node v) {
		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		// FIXME: bad accees to nodeToSuperNode
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			// add edge weight to supernode (self-loop) weight
			Gcon.setWeight(su, su, Gcon.weight(su) + G.weight(u, v));
		} else {
			// add edge weight to weight between two supernodes
			if (Gcon.hasEdge(su, sv)) {
				Gcon.setWeight(su, sv, Gcon.weight(su, sv) + G.weight(u, v));
			} else {
				// create new edge
				Gcon.insertEdge(su, sv, G.weight(u, v));
			}

		}
	}); // TODO: parallel?

	return std::make_pair(Gcon, nodeToSuperNode);

}

}
