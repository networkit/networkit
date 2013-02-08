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

	IndexMap<cluster, node> clusterToSuperNode(zeta.upperBound(), -1); // there is one supernode for each cluster

	// populate map cluster -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		if (! clusterToSuperNode.hasBeenSet(c)) {
			clusterToSuperNode[c] = Gcon.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});


	int64_t n = G.numberOfNodes();
	NodeMap<node> nodeToSuperNode(n);

	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = clusterToSuperNode[zeta.clusterOf(v)];
	});


	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.forWeightedEdges([&](node u, node v, edgeweight ew) {
		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		// FIXME: bad access to nodeToSuperNode
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			// add edge weight to supernode (self-loop) weight
			Gcon.setWeight(su, su, Gcon.weight(su, su) + ew);
		} else {
			// add edge weight to weight between two supernodes (or insert edge)
			Gcon.setWeight(su, sv, Gcon.weight(su, sv) + ew);
		}
	}); // TODO: parallel?

	return std::make_pair(Gcon, nodeToSuperNode);

}

}
