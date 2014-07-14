/*
 * ParallelPartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#include "ParallelPartitionCoarsening.h"
#include <omp.h>
#include "../auxiliary/Timer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {


std::pair<Graph, std::vector<node> > NetworKit::ParallelPartitionCoarsening::run(const Graph& G, const Partition& zeta) {

	Aux::Timer timer;
	timer.start();

	std::vector<node> subsetToSuperNode(zeta.upperBound(), none); // there is one supernode for each subset

	// populate map subset -> supernode
	node nextNodeId = 0;
	G.forNodes([&](node v){
		index c = zeta.subsetOf(v);
		if (subsetToSuperNode[c] == none) {
			subsetToSuperNode[c] = nextNodeId++;
		}
	});
	Graph Ginit(nextNodeId, true); // initial graph containing supernodes

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
	// combine local graphs in parallel
	Graph Gcombined(Ginit.numberOfNodes(), true, true); //
	Gcombined.parallelForNodes([&](node u) {
		for (index t = 0; t < nThreads; ++t) {
			localGraphs.at(t).forEdgesOf(u, [&](node u, node v) {
				Gcombined.addEdge(u, v, G.weight(u, v));
			});
		}
	});


	timer2.stop();
	INFO("combining coarse graphs took ", timer2.elapsedTag());

	timer.stop();
	INFO("parallel coarsening took ", timer.elapsedTag());

	return std::make_pair(localGraphs.back(), nodeToSuperNode);

}

} /* namespace NetworKit */
