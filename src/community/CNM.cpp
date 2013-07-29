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
	WARN("CNM code with PQ does not work properly yet!");


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
	std::vector<std::pair<double, index> > scores;
	std::map<index, Edge> indexToEdge; //!< necessary because PQ allows only integer values
	std::map<Edge, index> edgeToIndex;
	index edgeIndex = 0;
	G.forEdges([&](node u, node v) {
		double modScore = modScoring.edgeScore(u, v);
		scores.push_back(std::make_pair(-modScore, edgeIndex)); // negative due to minPQ
		indexToEdge.insert(std::make_pair(edgeIndex, std::make_pair(u, v)));
		edgeToIndex.insert(std::make_pair(std::make_pair(u, v), edgeIndex));
		++edgeIndex;
	});
	Aux::PriorityQueue<double, index> pq(scores);

	DEBUG("Initial PQ built before CNM loop");

	for (count clusters = n; clusters > 1; clusters--) {
		// determine best edge
		std::pair<double, index> pqMin = pq.extractMin();
		index bestEdgeIndex = pqMin.second;
		std::pair<node, node> bestEdge = indexToEdge[bestEdgeIndex];
		node best_u = bestEdge.first;
		node best_v = bestEdge.second;

		DEBUG("best edge score in iter " << (n - clusters) << ": " << - pqMin.first << ", edge: (" << best_u << "," << best_v << ")");

//		// remove dangling edges (those to be rewired) from PQ, then merge best edge
//		G.forEdgesOf(best_u, [&](node best_u, node neighbor) {
//			std::map<Edge, index>::iterator iter = edgeToIndex.find(std::make_pair(best_u, neighbor));
//			if (iter != edgeToIndex.end()) {
//				index idx = iter->second;
//				TRACE("removing index " << idx << ", edge: (" << best_u << "," << neighbor << ")");
//				std::pair<double, index> elem = std::make_pair(0.0, idx);
//				pq.remove(elem);
//			}
//		});
//
//		// remove dangling edges (those to be rewired) from PQ, then merge best edge
//		G.forEdgesOf(best_v, [&](node best_v, node neighbor) {
//			std::map<Edge, index>::iterator iter = edgeToIndex.find(std::make_pair(best_u, neighbor));
//			if (iter != edgeToIndex.end()) {
//				index idx = edgeToIndex[std::make_pair(best_v, neighbor)];
//				TRACE("removing index " << idx << ", edge: (" << best_v << "," << neighbor << ")");
//				std::pair<double, index> elem = std::make_pair(0.0, idx);
//				pq.remove(elem);
//			}
//		});

//		node newNode = G.mergeEdge(best_u, best_v, false);
//		TRACE("new node: " << newNode);


#if 0

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
#endif

		DEBUG("rewiring of edges incident to best_v and removing them from PQ");

		// rewiring of edges incident to best_v and removing them from PQ
		G.forWeightedEdgesOf(best_v, [&](node best_v, node neighbor, edgeweight ew) {
			pq.remove(std::make_pair(0.0, edgeToIndex[std::make_pair(best_v, neighbor)]));
			node endPoint = (best_v == neighbor) ? (best_u) : (neighbor);

			if (G.hasEdge(best_u, endPoint)) {
				// increase weight of existing edge
				G.setWeight(best_u, endPoint, ew + G.weight(best_u, endPoint));
			}
			else {
				// insert new edge to best_u
				G.addEdge(best_u, endPoint, ew);
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

		// update delta mod scores for incident edges and insert them into the PQ
		G.forEdgesOf(best_u, [&](node best_u, node neighbor) {
			if (best_u != neighbor) {
				double deltaMod = modScoring.edgeScore(best_u, neighbor);
				std::pair<double, index> elem = std::make_pair(-deltaMod, edgeIndex);
				indexToEdge.insert(std::make_pair(edgeIndex, std::make_pair(best_u, neighbor)));
				edgeToIndex.insert(std::make_pair(std::make_pair(best_u, neighbor), edgeIndex));
				TRACE("insert into PQ: (" << best_u << ", " << neighbor << ") with index " << edgeIndex << " and score " << deltaMod);
				++edgeIndex;
				pq.insert(elem);
			}
		});

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

