/*
 * Closeness.cpp
 *
 *  Created on: 03.10.2014
 *      Author: nemes
 */

#include <memory>
#include <queue>
#include <stack>

#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include "../components/ConnectedComponents.h"
#include "../components/StronglyConnectedComponents.h"
#include "../distance/BFS.h"
#include "../distance/Dijkstra.h"
#include "../distance/SSSP.h"
#include "Closeness.h"

namespace NetworKit {

Closeness::Closeness(const Graph &G, bool normalized,
                     const ClosenessVariant variant)
    : Centrality(G, normalized), n(G.upperNodeIdBound()), variant(variant) {
	if (variant == ClosenessVariant::standard) {
		checkConnectedComponents();
	}
}

Closeness::Closeness(const Graph &G, bool normalized, bool checkConnectedness)
    : Centrality(G, normalized), n(G.upperNodeIdBound()),
      variant(ClosenessVariant::standard) {
	if (checkConnectedness) {
		checkConnectedComponents();
	}
}

void Closeness::checkConnectedComponents() const {
	bool multipleComponents = false;
	if (G.isDirected()) {
		StronglyConnectedComponents scc(G);
		scc.run();
		multipleComponents = scc.numberOfComponents() > 1;
	} else {
		ConnectedComponents cc(G);
		cc.run();
		multipleComponents = cc.numberOfComponents() > 1;
	}
	if (multipleComponents) {
		throw std::runtime_error(
		    "Error: the standard definition of closeness is not defined on "
		    "disconnected graphs. On disconnected graphs, use the generalized "
		    "definition instead.");
	}
}

void Closeness::run() {
	scoreData.clear();
	scoreData.resize(n);
	if (variant == ClosenessVariant::generalized) {
		reachableNodes.assign(n, 0);
	}
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
		if (variant == ClosenessVariant::generalized) {
			reachableNodes[s] =
			    std::count_if(distances.begin(), distances.end(),
			                  [&](edgeweight dist) { return dist != infDist; });
		}

		double sum = 0;
		for (auto dist : distances) {
			if (dist != infDist) {
				sum += dist;
			}
		}

		scoreData[s] = sum == 0 ? 0
		                        : variant == ClosenessVariant::standard
		                              ? 1.0 / sum
		                              : (reachableNodes[s] - 1) / sum / (n - 1);
	});
	if (normalized) {
		G.parallelForNodes([&](node s) {
			scoreData[s] *=
			    (variant == ClosenessVariant::standard ? n : reachableNodes[s]) - 1;
		});
	}

	hasRun = true;
}

double Closeness::maximum() {
	return normalized ? 1.0f : (1.0f / (G.numberOfNodes() - 1));
}

} /* namespace NetworKit */
