/*
 * PartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#include "PartitionCoarsening.h"
#include <omp.h>


namespace NetworKit {


std::pair<Graph, std::vector<node> > NetworKit::PartitionCoarsening::run(Graph& G, Clustering& zeta) {

	Graph Gcon(0); // initial graph containing supernodes
	std::vector<node> subsetToSuperNode(zeta.upperBound(), none); // there is one supernode for each cluster

	// populate map subset -> supernode
	G.forNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		if (subsetToSuperNode[c] == none) {
			subsetToSuperNode[c] = Gcon.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});

	index z = G.upperNodeIdBound();
	std::vector<node> nodeToSuperNode(z);

	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = subsetToSuperNode[zeta.clusterOf(v)];
	});

	// make copies of initial graph
	int nThreads = omp_get_max_threads();
	std::vector<Graph> localGraphs(nThreads, Gcon);

	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.parallelForWeightedEdges([&](node u, node v, edgeweight ew) {
		int t = omp_get_thread_num();

		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];

		if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
			// add edge weight to supernode (self-loop) weight
			// TODO: Gcon.setWeight(su, su, Gcon.weight(su, su) + ew);
		} else {
			// add edge weight to weight between two supernodes (or insert edge)
		}
	});

	return std::make_pair(Gcon, nodeToSuperNode);

}

} /* namespace NetworKit */
