/*
 * Betweenness.cpp
 *
 *  Created on: 29.07.2014
 *      Author: cls, ebergamini
 */

#include <stack>
#include <queue>
#include <memory>
#include <omp.h>


#include "Betweenness.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/SignalHandling.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

Betweenness::Betweenness(const Graph& G, bool normalized, bool computeEdgeCentrality) : Centrality(G, normalized, computeEdgeCentrality) {

}

void Betweenness::run() {
	Aux::SignalHandler handler;
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);
	if (computeEdgeCentrality) {
		count z2 = G.upperEdgeIdBound();
		edgeScoreData.clear();
		edgeScoreData.resize(z2);
	}

	// thread-local scores for efficient parallelism
	count maxThreads = omp_get_max_threads();
	std::vector<std::vector<double> > scorePerThread(maxThreads, std::vector<double>(G.upperNodeIdBound()));
	DEBUG("score per thread: ", scorePerThread.size());
	DEBUG("G.upperEdgeIdBound(): ", G.upperEdgeIdBound());
	std::vector<std::vector<double> > edgeScorePerThread;
	if (computeEdgeCentrality) {
		edgeScorePerThread.resize(maxThreads, std::vector<double>(G.upperEdgeIdBound()));
	}
	DEBUG("edge score per thread: ", edgeScorePerThread.size());

	auto computeDependencies = [&](node s) {

		std::vector<double> dependency(z, 0.0);

		// run SSSP algorithm and keep track of everything
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, s, true, true));
		} else {
			sssp.reset(new BFS(G, s, true, true));
		}
		if (!handler.isRunning()) return;
		sssp->run();
		if (!handler.isRunning()) return;
		// compute dependencies for nodes in order of decreasing distance from s
		std::vector<node> stack = sssp->getStack();
		while (!stack.empty()) {
			node t = stack.back();
			stack.pop_back();
			for (node p : sssp->getPredecessors(t)) {
				// workaround for integer overflow in large graphs
				bigfloat tmp = sssp->numberOfPaths(p) / sssp->numberOfPaths(t);
				double weight;
				tmp.ToDouble(weight);
				double c= weight * (1 + dependency[t]);
				dependency[p] += c;
				if (computeEdgeCentrality) {
					edgeScorePerThread[omp_get_thread_num()][G.edgeId(p,t)] += c;
				}


			}
			if (t != s) {
				scorePerThread[omp_get_thread_num()][t] += dependency[t];
			}
		}
	};
	handler.assureRunning();
	G.balancedParallelForNodes(computeDependencies);
	handler.assureRunning();
	DEBUG("adding thread-local scores");
	// add up all thread-local values
	for (const auto &local : scorePerThread) {
		G.parallelForNodes([&](node v){
			scoreData[v] += local[v];
		});
	}
	if (computeEdgeCentrality) {
		for (const auto &local : edgeScorePerThread) {
			for (count i = 0; i < local.size(); i++) {
				edgeScoreData[i] += local[i];
			}
		}
	}
	if (normalized) {
		// divide by the number of possible pairs
		count n = G.numberOfNodes();
		count pairs = (n-2) * (n-1);
		count edges =  n    * (n-1);
		if (!G.isDirected()) {
			pairs = pairs / 2;
			edges = edges / 2;
		}
		G.forNodes([&](node u){
			scoreData[u] = scoreData[u] / pairs;
		});
		if (computeEdgeCentrality) {
			for (count edge = 0; edge < edgeScoreData.size(); edge++) {
				edgeScoreData[edge] =  edgeScoreData[edge] / edges;
			}
		}
	}

	hasRun = true;
}

double Betweenness::maximum(){
	if (normalized) {
		return 1;
	}
	double score;
	count n = G.numberOfNodes();
	if (G.isDirected()) {
		score = (n-1)*(n-2);
	} else {
		score = (n-1)*(n-2)/2;
	}
	return score;
}

} /* namespace NetworKit */
