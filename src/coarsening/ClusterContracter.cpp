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

Graph ClusterContracter::run(Graph& G, Clustering& zeta) {

	Graph* Gcon = new Graph();

	IndexMap<cluster, node> clusterToSuperNode(zeta.numberOfClusters()); // there is one supernode for each cluster

	//TODO: populate clusterToSuperNode
	G.forallNodes([&](node v){

	});

	// find supernode for node
	auto getSuperNode = [&](node v) {
		cluster cv = zeta[v];
		node sv = clusterToSuperNode[cv];
		return sv;
	};

	G.forallEdges([&](node u, node v) {
		node su = getSuperNode(u);
		node sv = getSuperNode(v);
		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			Gcon->setWeight(su, Gcon->weight(su) + G.weight(u, v));
		} else {
			Gcon->setWeight(su, sv, Gcon->weight(su, sv) + G.weight(u, v));
		}
	});

	return (*Gcon);
}

}
