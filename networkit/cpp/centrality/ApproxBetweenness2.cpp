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
#include "../auxiliary/SignalHandling.h"


#include <memory>

namespace NetworKit {

ApproxBetweenness2::ApproxBetweenness2(const Graph& G, count nSamples, bool normalized) : Centrality(G, normalized), nSamples(nSamples) {
}

void ApproxBetweenness2::run() {
	Aux::SignalHandler handler;

	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

	//std::vector<node> sampledNodes = G.nodes();
	std::vector<node> sampledNodes;

	// sample nodes
	for (count i = 0; i <= nSamples; ++i) {
	 	sampledNodes.push_back(G.randomNode());
	}

	for (node s : sampledNodes) {
		handler.assureRunning();
		// run single-source shortest path algorithm
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, s));
		} else {
			sssp.reset(new BFS(G, s));
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
			if (t == s){
				continue;
			}
			for (node p : sssp->getPredecessors(t)) {
				// TODO: make weighting factor configurable

				// workaround for integer overflow in large graphs
				bigfloat tmp = sssp->numberOfPaths(p) / sssp->numberOfPaths(t);
				double weight;
				tmp.ToDouble(weight);

				dependency[p] += (double(sssp->distance(p)) / sssp->distance(t)) * weight * (1 + dependency[t]);
			}
			scoreData[t] += dependency[t];
		}

	} // end for sampled nodes

	// extrapolate
	G.forNodes([&](node v) {
		scoreData[v] = scoreData[v] * (2 * G.numberOfNodes() / double(nSamples));
	});

	if (normalized) {
		// divide by the number of possible pairs
		count n = G.numberOfNodes();
		count pairs = (n-2) * (n-1);
		G.parallelForNodes([&](node u){
			scoreData[u] = scoreData[u] / pairs;
		});
	}

	hasRun = true;
}


} /* namespace NetworKit */
