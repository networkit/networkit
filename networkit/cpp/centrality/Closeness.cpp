/*
 * Closeness.cpp
 *
 *  Created on: 03.10.2014
 *      Author: nemes
 */

#include <stack>
#include <queue>
#include <memory>

#include "Closeness.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

Closeness::Closeness(const Graph& G, bool normalized) : Centrality(G, normalized, false) {

}

void Closeness::run() {
	// TODO: we might want to add a parallel version
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

	G.forNodes([&](node s) {
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, s, true, true));
		} else {
			sssp.reset(new BFS(G, s, true, true));
		}
		sssp->run();

		std::vector<edgeweight> distances = sssp->getDistances();
		double sum = 0;
		for (auto dist : distances) {
			sum += dist;
		}
		scoreData[s] = 1 / sum;

	});
	if (normalized) {
		//TODO
	}
}


} /* namespace NetworKit */
