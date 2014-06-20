/*
* ApproxBetweenness2.cpp
*
*  Created on: 13.06.2014
*      Author: Christian Staudt, Elisabetta Bergamini
*/


#include "ApproxBetweenness2.h"
#include "../graph/BFS.h"
#include "../graph/Dijkstra.h"
#include "../graph/SSSP.h"

namespace NetworKit {

ApproxBetweenness2::ApproxBetweenness2(const Graph& G, count nSamples, bool normalized) : Centrality(G, normalized), nSamples(nSamples) {
}

void ApproxBetweenness2::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

	std::vector<node> sampledNodes = G.nodes();

	// // sample nodes
	// for (count i = 0; i <= nSamples; ++i) {
	// 	sampledNodes.push_back(G.randomNode());
	// }

	for (node s : sampledNodes) {
		// run single-source shortest path algorithm
		SSSP* sssp;
		if (G.isWeighted()) {
			sssp = new Dijkstra(G, s);
		} else {
			sssp = new BFS(G, s);
		}

		sssp->run();

		// create stack of nodes in non-decreasing order of distance
		std::vector<node> stack = G.nodes();
		std::sort(stack.begin(), stack.end(), [&](node u, node v){
			return (sssp->distance(u) > sssp->distance(v));
		});

		// compute dependencies and add the contributions to the centrality score
		std::vector<double> dependency(G.upperNodeIdBound(), 0.0);
		for (node t : stack) {
			for (node p : sssp->getPredecessors(t)) {
				// TODO: make weighting factor configurable
				// (sssp->distance(p) / sssp->distance(t))
				dependency[p] += (double(sssp->numberOfPaths(p)) / sssp->numberOfPaths(t)) * (1 + dependency[t]);
			}

			scoreData[t] += dependency[t];
		}

	} // end for sampled nodes

	// extrapolate
	//G.forNodes([&](node v) {
	//	scoreData[v] = scoreData[v] * (G.numberOfNodes() / double(nSamples));
	//});

	if (normalized) {
		// divide by the number of possible pairs
		count n = G.numberOfNodes();
		count pairs = (n-2) * (n-1);
		G.parallelForNodes([&](node u){
			scoreData[u] = scoreData[u] / pairs;
		});
	}

}


} /* namespace NetworKit */
