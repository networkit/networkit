/*
 * ParallelPartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#include "ParallelPartitionCoarsening.h"
#include <omp.h>
#include "../graph/GraphBuilder.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ParallelPartitionCoarsening::ParallelPartitionCoarsening(const Graph& G, const Partition& zeta, bool useGraphBuilder) : GraphCoarsening(G), zeta(zeta),	useGraphBuilder(useGraphBuilder) {

}

void ParallelPartitionCoarsening::run() {
	Aux::Timer timer;
	timer.start();

	Partition nodeToSuperNode = zeta;
	nodeToSuperNode.compact((zeta.upperBound() <= G.upperNodeIdBound())); // use turbo if the upper id bound is <= number of nodes
	count nextNodeId = nodeToSuperNode.upperBound();

	Graph Gcombined;
	if (!useGraphBuilder) {
		Graph Ginit(nextNodeId, true); // initial graph containing supernodes

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
		// Graph Gcombined(Ginit.numberOfNodes(), true); //
		Gcombined = Graph(Ginit.numberOfNodes(), true);

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

		assert (Gcombined.checkConsistency());

		// stop both timers before printing
		timer2.stop();
		INFO("combining coarse graphs took ", timer2.elapsedTag());
	} else {
		std::vector< std::vector<node> > nodesPerSuperNode(nextNodeId);
		G.forNodes([&](node v) {
			node sv = nodeToSuperNode[v];
			nodesPerSuperNode[sv].push_back(v);
		});

		// iterate over edges of G and create edges in coarse graph or update edge and node weights in Gcon
		DEBUG("create edges in coarse graphs");
		GraphBuilder b(nextNodeId, true, false);
		#pragma omp parallel for schedule(guided)
		for (node su = 0; su < nextNodeId; su++) {
			std::map<index, edgeweight> outEdges;
			for (node u : nodesPerSuperNode[su]) {
				G.forNeighborsOf(u, [&](node v, edgeweight ew) {
					node sv = nodeToSuperNode[v];
					if (su != sv || u >= v) { // count edges inside uv only once (we iterate over them twice)
						outEdges[sv] += ew;
					}
				});
			}
			for (auto it : outEdges) {
				b.addHalfEdge(su, it.first, it.second);
			}
		}

		Gcombined = b.toGraph(false);
	}

	timer.stop();
	INFO("parallel coarsening took ", timer.elapsedTag());
	Gcoarsened = std::move(Gcombined);
	nodeMapping = nodeToSuperNode.getVector();
	hasRun = true;
}

} /* namespace NetworKit */
