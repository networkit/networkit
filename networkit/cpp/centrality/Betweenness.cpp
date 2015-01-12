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
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

Betweenness::Betweenness(const Graph& G, bool normalized) : Centrality(G, normalized) {

}

void Betweenness::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

	// thread-local scores for efficient parallelism
	count maxThreads = omp_get_max_threads();
	std::vector<std::vector<double> > scorePerThread(maxThreads, std::vector<double>(G.upperNodeIdBound()));


	auto computeDependencies = [&](node s) {

		std::vector<double> dependency(z, 0.0);

		// run SSSP algorithm and keep track of everything
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, s, true, true));
		} else {
			sssp.reset(new BFS(G, s, true, true));
		}

		sssp->run();

		// compute dependencies for nodes in order of decreasing distance from s
		std::stack<node> stack = sssp->getStack();
		while (!stack.empty()) {
			node t = stack.top();
			stack.pop();
			for (node p : sssp->getPredecessors(t)) {
				// workaround for integer overflow in large graphs
				bigfloat tmp = sssp->numberOfPaths(p) / sssp->numberOfPaths(t);
				double weight;
				tmp.ToDouble(weight);

				dependency[p] += weight * (1 + dependency[t]);
			}
			if (t != s) {
				scorePerThread[omp_get_thread_num()][t] += dependency[t];
			}
		}
	};



	G.balancedParallelForNodes(computeDependencies);

	INFO("adding thread-local scores");
	// add up all thread-local values
	for (auto local : scorePerThread) {
		G.parallelForNodes([&](node v){
			scoreData[v] += local[v];
		});
	}

	if (normalized) {
		// divide by the number of possible pairs
		count n = G.numberOfNodes();
		count pairs = (n-2) * (n-1);
		G.forNodes([&](node u){
			scoreData[u] = scoreData[u] / pairs;
		});
	}
}






} /* namespace NetworKit */
