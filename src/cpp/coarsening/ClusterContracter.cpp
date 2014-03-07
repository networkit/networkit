/*
 * ClusterContractor.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusterContractor.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

ClusterContractor::ClusterContractor() {

}

ClusterContractor::~ClusterContractor() {

}

std::pair<Graph, std::vector<node> > ClusterContractor::run(const Graph& G, const Partition& zeta) {

	Aux::Timer timer;
	timer.start();

	Graph Gcon(0, true); // empty weighted graph

	std::vector<node> clusterToSuperNode(zeta.upperBound(), none); // there is one supernode for each cluster

	DEBUG("populate map cluster -> supernode");
	G.forNodes([&](node v){
		index c = zeta.subsetOf(v);//zeta.clusterOf(v);Cluster
		if (clusterToSuperNode[c] == none) {
			clusterToSuperNode[c] = Gcon.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});

	index z = G.upperNodeIdBound();
	std::vector<node> nodeToSuperNode(z, none);

	DEBUG("set entries node -> supernode");
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = clusterToSuperNode[zeta.subsetOf(v)];
	});


	DEBUG("iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon");
	G.forWeightedEdges([&](node u, node v, edgeweight ew) {
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
