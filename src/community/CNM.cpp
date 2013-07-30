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
	typedef std::pair<node, node> Edge;

	// copy graph because we make changes due to merges
	Graph G = graph;
	count n = G.numberOfNodes();

	// start with singleton clustering
	Clustering clustering(n);
	clustering.allToSingletons();

	// record current modularity
	Modularity modularityInspector;
	double bestModularity = modularityInspector.getQuality(clustering, G);
	Clustering bestClustering = clustering;

	// insert edges into priority queue, key: delta mod score
	ModularityScoring<double> modScoring(G);
	std::vector<std::pair<double, std::pair<node, node> > > scores;
//	std::map<index, Edge> indexToEdge; //!< necessary because PQ allows only integer values
//	std::map<Edge, index> edgeToIndex;
//	index edgeIndex = 0;
	G.forEdges([&](node u, node v) {
		double modScore = modScoring.edgeScore(u, v);
		std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
		scores.push_back(std::make_pair(-modScore, edge)); // negative due to minPQ
//		indexToEdge.insert(std::make_pair(edgeIndex, std::make_pair(u, v)));
//		edgeToIndex.insert(std::make_pair(std::make_pair(u, v), edgeIndex));
//		++edgeIndex;
	});

	Aux::PriorityQueue<double, std::pair<node, node> > pq(scores);

	DEBUG("Initial PQ built before CNM loop");

	while (pq.size() > 0) {

		// determine best edge
		std::pair<double, std::pair<node, node> > pqMin = pq.extractMin();
		node best_u = (pqMin.second).first;
		node best_v = (pqMin.second).second;

		DEBUG("best_u: " << best_u << ", best_v: " << best_v);
		DEBUG("rewiring of edges incident to best_v and removing them from PQ");

		// remove best_u edges from PQ (has to be done first)
		G.forWeightedEdgesOf(best_u, [&](node best_u, node neighbor, edgeweight ew) {
			std::pair<node, node> edge = (best_u <= neighbor) ? (std::make_pair(best_u, neighbor)) : (std::make_pair(neighbor, best_u));
			pq.remove(std::make_pair(0.0, edge));
		});

		// rewiring of edges incident to best_v and removing them from PQ
		G.forWeightedEdgesOf(best_v, [&](node best_v, node neighbor, edgeweight ew) {
			std::pair<node, node> edge = (best_v <= neighbor) ? (std::make_pair(best_v, neighbor)) : (std::make_pair(neighbor, best_v));
			node endPoint = (best_v == neighbor) ? (best_u) : (neighbor);
			if (best_v != neighbor) {
				pq.remove(std::make_pair(0.0, edge));
			}

			if (G.hasEdge(best_u, endPoint)) {
				// increase weight of existing edge
				G.increaseWeight(best_u, endPoint, ew);
				TRACE("increase weight of edge " << best_u << "," << endPoint << " by " << ew);
			}
			else {
				// insert new edge to best_u
				G.addEdge(best_u, endPoint, ew);
				TRACE("insert edge " << best_u << "," << endPoint);
			}
		});

		DEBUG("remove best_v from graph");

		// remove best_v from graph
		G.forEdgesOf(best_v, [&](node best_v, node neighbor) {
			G.removeEdge(best_v, neighbor);
		});
		G.removeNode(best_v);

		// merge clusters best_u and best_v
		clustering.mergeClusters(clustering.clusterOf(best_u), clustering.clusterOf(best_v));

#if 0
		G.forEdges([&](node u, node v) {
			std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
			pq.remove(std::make_pair(0.0, edge));
			double deltaMod = modScoring.edgeScore(best_u, neighbor);
			pq.insert(std::make_pair(deltaMod, edge));
		});
#else
		// update delta mod scores for incident edges and insert them into the PQ, needs to
		// done after all edges have been rewired
		G.forEdgesOf(best_u, [&](node best_u, node neighbor) {
			if (best_u != neighbor) {
				std::pair<node, node> edge = (best_u <= neighbor) ? (std::make_pair(best_u, neighbor)) : (std::make_pair(neighbor, best_u));
				double deltaMod = modScoring.edgeScore(best_u, neighbor);
				TRACE("insert into PQ: " << -deltaMod << " for " << best_u << " and " << neighbor);
				std::pair<double, std::pair<node, node> > elem = std::make_pair(-deltaMod, edge);
				pq.insert(elem);
			}
		});
#endif

		// compute new modularity
		double newModularity = modularityInspector.getQuality(clustering, graph);

		// record best solutions
		if (newModularity > bestModularity) {
			TRACE("new mod value: " << newModularity);
			bestModularity = newModularity;
			bestClustering = clustering;
		}
	}

	return bestClustering;
}

} // namespace

