/*
 * CNM.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: Matthias Wolf & Michael Wegner
 */

#include "CNM.h"
#include "../scoring/ModularityScoring.h"
#include "../clustering/Modularity.h"
#include <limits>

namespace NetworKit{

CNM::CNM() {

}

CNM::~CNM() {

}

Clustering CNM::run(Graph &graph) {
	// copy graph because we make changes due to merges
	Graph G = graph;
	count n = G.numberOfNodes();

	// start with singleton clustering
	Clustering clustering(n);
	clustering.allToSingletons();

	// record current modularity
	Modularity modularityInspector;
	double bestModularity = modularityInspector.getQuality(clustering, G);
	Clustering bestClustering(n);

#if 0

	// insert edges into priority queue, key: delta mod score
	ModularityScoring<double> modScoring(G);
	std::vector<double, std::pair<node, node> > scores;
	G.forEdges([&](node u, node v) {
		scores.push_back(modScoring.edgeScore(u, v), std::make_pair(u, v));
	});
	Aux::PriorityQueue<double, index> pq(scores);

	for (count clusters = n; clusters > 1; clusters--) {
		// determine best edge
		std::pair<node, node> bestEdge = pq.extractMin().second;
		node best_u = bestEdge.first;
		node best_v = bestEdge.second;

		// TODO (in graph): merge edge to supernode and adjust neighborhood
		// mergeEdge(node u, node v, bool discardSelfLoop = false)



		double bestDelta = std::numeric_limits<double>::lowest();
		node best_u, best_v;

		// find edge to be merged (highest modularity score)
		auto calcBestDelta = [&](node u, node v) {
			if (u == v) return;

			double delta = modScoring.edgeScore(u, v);
			if (delta > bestDelta) {
				bestDelta = delta;
				best_u = u;
				best_v = v;
			}
		};

		// compute modularity scores for all edges
		G.forEdges(calcBestDelta);

		// rewiring of edges
		auto changeMountPoint = [&G, &best_u, &best_v](node u, node v, edgeweight ew) {
//			G.removeEdge(u, v);
			assert(best_v == u);

			node endPoint = v;
			if (u == v) {
				// loop
				endPoint = best_u;
			}

			if (G.hasEdge(best_u, endPoint)) {
				G.setWeight(best_u, endPoint, ew + G.weight(best_u, endPoint));
			} else {
				G.addEdge(best_u, endPoint, ew);
			}
		};

		// rewiring
		G.forWeightedEdgesOf(best_v, changeMountPoint);
		auto removeEdges = [&G](node u, node v) {
			G.removeEdge(u, v);
		};

		// best_u remains => remove dangling edges incident to best_v
		G.forEdgesOf(best_v, removeEdges);

		// merge clusters best_u and best_v
		clustering.mergeClusters(clustering.clusterOf(best_u), clustering.clusterOf(best_v));

		// compute new modularity
		double newModularity = modularityInspector.getQuality(clustering, G);

		// record best solutions
		if (newModularity > bestModularity) {
			TRACE("new mod value: " << newModularity);
			bestModularity = newModularity;
			bestClustering = clustering;
		}
	}

#endif

	return bestClustering;
}

} // namespace

