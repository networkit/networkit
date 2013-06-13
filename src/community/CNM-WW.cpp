/*
 * CNM.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: Matthias Wolf & Michael Wegner
 */

#include "CNM-WW.h"
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

	for (count clusters = n; clusters > 1; clusters--) {
		ModularityScoring<double> modScoring(G);
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

		// ???
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

		// ???
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
			DEBUG(newModularity);
			bestModularity = newModularity;
			bestClustering = clustering;
		}
	}

	return bestClustering;
}

}

