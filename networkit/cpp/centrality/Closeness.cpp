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
#include "../components/ConnectedComponents.h"


namespace NetworKit {

Closeness::Closeness(const Graph& G, bool normalized, bool checkConnectedness) : Centrality(G, normalized) {
	// TODO: extend closeness definition to make check for connectedness unnecessary
	if (checkConnectedness) {
		ConnectedComponents compo(G);
		compo.run();
		if (compo.numberOfComponents() != 1) {
			throw std::runtime_error("Closeness is undefined on disconnected graphs");
		}
	}
}

void Closeness::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);
	edgeweight infDist = std::numeric_limits<edgeweight>::max();

	G.parallelForNodes([&](node s) {
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
			if (dist != infDist ) {
				sum += dist;
			}
		}
		scoreData[s] = 1 / sum;

	});
	if (normalized) {
		G.forNodes([&](node u){
			scoreData[u] = scoreData[u] * (G.numberOfNodes() - 1);
		});
	}

	hasRun = true;
}

double Closeness::maximum() {
	return (double) 1 / (G.numberOfNodes() - 1);
}

} /* namespace NetworKit */
