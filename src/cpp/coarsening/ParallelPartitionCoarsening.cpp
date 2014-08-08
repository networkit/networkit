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

	DEBUG("populate map subset -> supernode");
	node nextNodeId = 0;
	G.forNodes([&](node v){
		index c = zeta.subsetOf(v);
		if (subsetToSuperNode[c] == none) {
			subsetToSuperNode[c] = nextNodeId++;
		}
	});
	Graph Ginit(nextNodeId, true); // initial graph containing supernodes

	index z = G.upperNodeIdBound();
	std::vector<node> nodeToSuperNode(z, none);

	// set entries node -> supernode
	DEBUG("set entries node -> supernode");
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = subsetToSuperNode[zeta.subsetOf(v)];
	});

	// make copies of initial graph
	count nThreads = omp_get_max_threads();
	std::vector<Graph> localGraphs(nThreads, Ginit); // thread-local graphs


	// iterate over edges of G and create edges in coarse graph or update edge and node weights in Gcon
	DEBUG("create edges in coarse graphs");
	G.parallelForEdges([&](node u, node v, edgeweight ew) {
		index t = omp_get_thread_num();

		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		localGraphs.at(t).increaseWeight(su, sv, ew);

	});


	Aux::Timer timer2;
	timer2.start();
	// combine local graphs in parallel
	Graph Gcombined(Ginit.numberOfNodes(), true); //

	std::vector<count> numEdges(nThreads);


	// access internals of Graph to write adjacencies
	auto threadSafeIncreaseWeight = [&](node u, node v, edgeweight ew) {

		index vi = Gcombined.indexInOutEdgeArray(u, v);
		if (vi == none) {
			index t = omp_get_thread_num();
			if (u == v) {
				numEdges[t] += 2;
			} else {
				numEdges[t] += 1; // normal edges count half
			}
			Gcombined.outDeg[u]++;
			Gcombined.outEdges[u].push_back(v);
			Gcombined.outEdgeWeights[u].push_back(ew);
		} else {
			Gcombined.outEdgeWeights[u][vi] += ew;
		}

	};

	DEBUG("combining graphs");
	Gcombined.balancedParallelForNodes([&](node u) {
		for (index l = 0; l < nThreads; ++l) {
			localGraphs.at(l).forEdgesOf(u, [&](node u, node v, edgeweight w) {
				TRACE("increasing weight of (", u, v, ") to", w);
				threadSafeIncreaseWeight(u, v, w);
			});
		}
	});


	// ensure consistency of data structure
	DEBUG("numEdges: ", numEdges);
	count twiceM = std::accumulate(numEdges.begin(), numEdges.end(), 0);
	assert (twiceM % 2 == 0);
	Gcombined.m = (twiceM / 2);

	assert (G.consistencyCheck());

	// stop both timers before printing
	timer2.stop();
	timer.stop();
	INFO("combining coarse graphs took ", timer2.elapsedTag());
	INFO("parallel coarsening took ", timer.elapsedTag());

	return {std::move(Gcombined), std::move(nodeToSuperNode)};

}

} /* namespace NetworKit */
