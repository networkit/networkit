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
#include "../auxiliary/PrioQueue.h"

namespace NetworKit{

CNM::CNM() {

}

CNM::~CNM() {

}

Partition CNM::run(Graph &graph) {
	// copy graph because we make changes due to merges
	Graph G = graph;
	count n = G.numberOfNodes();

	// start with singleton clustering
	Partition clustering(n);
	clustering.allToSingletons();

	// record current modularity
	Modularity modularityInspector;
	double bestModularity = modularityInspector.getQuality(clustering, G);
	Partition bestClustering = clustering;

	// insert edges into priority queue, key: delta mod score
	ModularityScoring<double> modScoring(G);
	std::vector<double> scores;
	std::map<std::pair<node, node>, index> edgeToIndex;
	std::map<index, std::pair<node, node> > indexToEdge;
	index idx = 0;
	G.forEdges([&](node u, node v) {
		// get modularity score
		double modScore = modScoring.edgeScore(u, v);

		// create edge entry (unordered pair!)
		std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
		scores.push_back(-modScore); // negative due to minPQ (and we are maximizing here)

		// enable mapping between edge and edge index -- and vice versa
		edgeToIndex.insert(std::make_pair(edge, idx));
		indexToEdge.insert(std::make_pair(idx, edge));
		++idx;
	});
	Aux::PrioQueue<double, index> pq(scores);

	TRACE("Initial PQ built before CNM loop");

	count iterations = 0;
	while (iterations < n-1) {
		// determine best edge
		node best_u, best_v;
		do {
			std::pair<double, index> pqMin = pq.extractMin();
			std::pair<node, node> currentEdge = indexToEdge[pqMin.second];
			best_u = currentEdge.first;
			best_v = currentEdge.second;
		} while (! G.hasEdge(best_u, best_v));

		TRACE("best_u: " , best_u , ", best_v: " , best_v);
//		TRACE("rewiring of edges incident to best_v and removing them from PQ");

		// merge clusters best_u and best_v
		index newCluster = clustering.mergeSubsets(clustering.subsetOf(best_u), clustering.subsetOf(best_v));

//		TRACE("merged clustering, about to merge nodes");

		node newNode = G.mergeEdge(best_u, best_v, false);
		assert(newNode != none);
		clustering.extend();
		clustering[newNode] = newCluster;

//		TRACE("graph changed, about to update PQ");

		G.forEdgesOf(best_u, [&](node best_u, node neighbor) {
			std::pair<node, node> edge = (best_u <= neighbor) ? (std::make_pair(best_u, neighbor)) : (std::make_pair(neighbor, best_u));
			pq.remove(edgeToIndex[edge]);
		});
		G.forEdgesOf(best_v, [&](node best_v, node neighbor) {
			std::pair<node, node> edge = (best_v <= neighbor) ? (std::make_pair(best_v, neighbor)) : (std::make_pair(neighbor, best_v));
			pq.remove(edgeToIndex[edge]);
		});


		// update delta mod scores for incident edges and insert them into the PQ, needs to
		// done after all edges have been rewired
		G.forEdgesOf(newNode, [&](node newNode, node neighbor) {
			if (newNode != neighbor) {
				ModularityScoring<double> modScoring(G);
				double score = modScoring.edgeScore(newNode, neighbor);
//				TRACE("score of new edge: ", score);

				std::pair<node, node> edge = std::make_pair(newNode, neighbor);
				edgeToIndex.insert(std::make_pair(edge, idx));
				indexToEdge.insert(std::make_pair(idx, edge));
				++idx;

//				TRACE("insert (-score, idx): ", -score, ", ", idx);
				pq.insert(-score, idx);
			}
		});

//		TRACE("updated PQ, about to compute quality");

		// compute new modularity
		double newModularity = modularityInspector.getQuality(clustering, graph);
		TRACE("pq size: ", pq.size(), ", current mod value: " , newModularity);

		// record best solutions
		if (newModularity > bestModularity) {
			DEBUG("ha, improved so far best mod value: " , bestModularity);
			bestModularity = newModularity;
			bestClustering = clustering;
		}

		++iterations;
	}

	return bestClustering;
}

} // namespace

