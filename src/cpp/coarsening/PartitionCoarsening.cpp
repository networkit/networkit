/*
 * PartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#include "PartitionCoarsening.h"
#include <omp.h>
#include "../auxiliary/Timer.h"


namespace NetworKit {


std::pair<Graph, std::vector<node> > NetworKit::PartitionCoarsening::run(const Graph& G, const Partition& zeta) {

	Aux::Timer timer;
	timer.start();

	Graph Ginit(0, true); // initial graph containing supernodes
	std::vector<node> subsetToSuperNode(zeta.upperBound(), none); // there is one supernode for each cluster

	// populate map subset -> supernode
	G.forNodes([&](node v){
		index c = zeta.subsetOf(v);
		if (subsetToSuperNode[c] == none) {
			subsetToSuperNode[c] = Ginit.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});

	index z = G.upperNodeIdBound();
	std::vector<node> nodeToSuperNode(z);

	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = subsetToSuperNode[zeta.subsetOf(v)];
	});

	// make copies of initial graph
	count nThreads = omp_get_max_threads();
	std::vector<Graph> localGraphs(nThreads, Ginit); // thread-local graphs

	// iterate over edges of G and create edges in coarse graph or update edge and node weights in Gcon
	G.parallelForWeightedEdges([&](node u, node v, edgeweight ew) {
		index t = omp_get_thread_num();

		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		localGraphs.at(t).increaseWeight(su, sv, ew);

	});

	Aux::Timer timer2;
	timer2.start();
	// combine local graphs
	for (index i = 0; i < (nThreads - 1); ++i) {
		localGraphs[i].forWeightedEdges([&](node u, node v, edgeweight ew) {
			localGraphs.at(i+1).increaseWeight(u, v, ew);
		});
	}
	timer2.stop();
	INFO("combining coarse graphs took ", timer2.elapsedTag());

	timer.stop();
	INFO("parallel coarsening took ", timer.elapsedTag());

	return std::make_pair(localGraphs.back(), nodeToSuperNode);

}

} /* namespace NetworKit */
