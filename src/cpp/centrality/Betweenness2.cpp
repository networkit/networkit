/*
 * Betweenness.cpp
 *
 *  Created on: 19.02.2014
 *      Author: hm
 */

#include <stack>
#include <queue>
#include <memory>

#include "Betweenness2.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

Betweenness2::Betweenness2(const Graph& G, bool normalized) : Centrality(G, normalized) {

}

void Betweenness2::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

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
				dependency[p] += (double(sssp->numberOfPaths(p)) / sssp->numberOfPaths(t)) * (1 + dependency[t]);
			}
			if (t != s) {
				scoreData[t] += dependency[t];
			}
		}
	};

	G.forNodes(computeDependencies);
}






} /* namespace NetworKit */
