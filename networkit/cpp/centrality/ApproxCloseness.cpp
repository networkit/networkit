/*
* ApproxCloseness.cpp
*
*  Created on: 16.06.2015
*      Author: Arie Slobbe
*/


#include "ApproxCloseness.h"
#include "../graph/BFS.h"
#include "../graph/Dijkstra.h"
#include "../graph/SSSP.h"
#include "../auxiliary/SignalHandling.h"


#include <memory>

namespace NetworKit {

ApproxCloseness::ApproxCloseness(const Graph& G, count nSamples, bool normalized) : Centrality(G, normalized), nSamples(nSamples) {
	if (G.isDirected()) {
		throw std::runtime_error("Graph is directed. ApproxCloseness only runs on undirected graphs.");
	}
}

void ApproxCloseness::run() {
	Aux::SignalHandler handler;

	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

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

		// increment scoreData with SSSP values
		std::vector<edgeweight> distances = sssp->getDistances();
		DEBUG("distances: ", distances);
		G.forNodes([&](node u){
			scoreData[u] += distances[u];
		});
		DEBUG("scoreData[", s, "]: ", scoreData[s]);


	} // end for sampled nodes

	// Compute estimated closeness centrality scores.
  count nNodes = G.numberOfNodes();
	if (normalized) {
		G.parallelForNodes([&](node u){
			scoreData[u] = (nSamples * (nNodes - 1)) / (scoreData[u] * nNodes);
		});
	} else {
		G.parallelForNodes([&](node u){
			scoreData[u] = nSamples / (scoreData[u] * nNodes);
		});
	}

	hasRun = true;
}

double ApproxCloseness::maximum() {
	return (double) 1 / (G.numberOfNodes() - 1);
}

} /* namespace NetworKit */
