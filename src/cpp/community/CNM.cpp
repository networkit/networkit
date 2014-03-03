/*
 * CNM.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: Matthias Wolf & Michael Wegner
 */

#include "CNM.h"
#include "../scoring/ModularityScoring.h"
#include "../community/Modularity.h"
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

	auto createEdge([&](node u, node v) {
		std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
		return edge;
	});


	// start with singleton clustering
	Partition clustering(n);
	clustering.allToSingletons();

	// record current modularity
	double gTotalEdgeWeight = G.totalEdgeWeight();
	Modularity modularityInspector;
	modularityInspector.setTotalEdgeWeight(gTotalEdgeWeight);
	double bestModularity = modularityInspector.getQuality(clustering, G);
	Partition bestClustering = clustering;

	// insert edges into priority queue, key: delta mod score
	ModularityScoring<double> modScoring(G);
	std::vector<double> scores;
	std::map<std::pair<node, node>, index> edgeToIndex;
	std::map<index, std::pair<node, node> > indexToEdge;
	index idx = 0;
	G.forEdges([&](node u, node v) {
		if (u != v) {
			// get modularity score
			double modScore = modScoring.edgeScore(u, v);

			// create edge entry (unordered pair!)
			std::pair<node, node> edge = createEdge(u, v);
			scores.push_back(-modScore); // negative due to minPQ (and we are maximizing here)

			// enable mapping between edge and edge index -- and vice versa
			edgeToIndex.insert(std::make_pair(edge, idx));
			indexToEdge.insert(std::make_pair(idx, edge));
			++idx;
		}
	});
	Aux::PrioQueue<double, index> pq(scores);
//	TRACE("Initial PQ built before CNM loop");
//	pq.print();

	// free mem
	scores.clear();

	count iterations = 0;
	while ((iterations < n-1) && (pq.size() > 0)) {
		// determine best edge
		std::pair<double, index> pqOpt = pq.extractMin();
		std::pair<node, node> currentEdge = indexToEdge[pqOpt.second];
		node best_u = currentEdge.first;
		node best_v = currentEdge.second;
//		TRACE("best_u: " , best_u , ", best_v: " , best_v, ", ew: ", G.weight(best_u, best_v), ", deltaMod: ", -pqOpt.first);

		// merge clusters best_u and best_v
		index newCluster = clustering.mergeSubsets(clustering.subsetOf(best_u), clustering.subsetOf(best_v));

		// adapt PQ accordingly
		G.forEdgesOf(best_u, [&](node best_u, node neighbor) {
			pq.remove(edgeToIndex[createEdge(best_u, neighbor)]);
		});
		G.forEdgesOf(best_v, [&](node best_v, node neighbor) {
			pq.remove(edgeToIndex[createEdge(best_v, neighbor)]);
		});

		// merge incident nodes in G
		node newNode = G.mergeEdge(best_u, best_v, false);
		assert(newNode != none);

		// adapt clustering accordingly
		clustering.extend();
		clustering[newNode] = newCluster;

		// update delta mod scores for incident edges and insert them into the PQ,
		// needs to done after all edges have been rewired
		G.forEdgesOf(newNode, [&](node newNode, node neighbor) {
			if (newNode != neighbor) {
				// determine edge score
//				assert(G.totalEdgeWeight() == gTotalEdgeWeight);
				ModularityScoring<double> modScoring(G, gTotalEdgeWeight);
				double score = modScoring.edgeScore(newNode, neighbor);

				// insert edge score and corresponding ID into PQ
				std::pair<node, node> edge = createEdge(newNode, neighbor);
				edgeToIndex.insert(std::make_pair(edge, idx));
				indexToEdge.insert(std::make_pair(idx, edge));

//				TRACE("insert ", -score, ", ", idx, " into PQ for edge ", newNode, "-", neighbor);
				pq.insert(-score, idx);

				++idx;
			}
		});

		// compute new modularity
//		assert(G.totalEdgeWeight() == gTotalEdgeWeight);
		modularityInspector.setTotalEdgeWeight(gTotalEdgeWeight);
		double newModularity = modularityInspector.getQuality(clustering, graph);
//		TRACE("pq size: ", pq.size(), ", current mod value: " , newModularity);

		// record best solutions
		if (newModularity > bestModularity) {
//			DEBUG("ha, improved so far best mod value: " , bestModularity);
			bestModularity = newModularity;
			bestClustering = clustering;
		}

		++iterations;
//		pq.print();
	}

	return bestClustering;
}

} // namespace

