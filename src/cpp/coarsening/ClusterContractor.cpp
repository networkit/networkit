/*
 * ClusterContractor.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusterContractor.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ClusterContractor::ClusterContractor() {

}

ClusterContractor::~ClusterContractor() {

}

std::pair<Graph, std::vector<node> > ClusterContractor::run(const Graph& G, const Partition& zeta) {

	Aux::Timer timer;
	timer.start();


	std::vector<node> clusterToSuperNode(zeta.upperBound()+1, none); // there is one supernode for each cluster
	// +1 is an experimental fix

	node nextNodeId = 0;
	DEBUG("populate map cluster -> supernode");
	G.forNodes([&](node v){
		index c = zeta.subsetOf(v);
		if (clusterToSuperNode[c] == none) {
			clusterToSuperNode[c] = nextNodeId++;
		}
	});
	Graph Gcon(nextNodeId, true); // empty weighted graph

	index z = G.upperNodeIdBound() + 1; // +1 is an experimental fix
	std::vector<node> nodeToSuperNode(z, none);

	DEBUG("set entries node -> supernode");
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = clusterToSuperNode[zeta.subsetOf(v)];
	});


	DEBUG("iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon");
	G.forEdges([&](node u, node v, edgeweight ew) {
		// TRACE(Gcon.upperNodeIdBound()," ",nodeToSuperNode[u]," ",nodeToSuperNode[v]);
		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		//TRACE("edge (", su, ", ", sv, ")");
		// add edge weight to weight between two supernodes (or insert edge)
		Gcon.increaseWeight(su, sv, ew);
	});

	timer.stop();
	INFO("sequential coarsening took ", timer.elapsedTag());

	return std::make_pair(Gcon, nodeToSuperNode);

}

}
