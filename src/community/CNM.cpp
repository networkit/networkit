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
	Clustering bestClustering = clustering;

	// insert edges into priority queue, key: delta mod score
	ModularityScoring<double> modScoring(G);
	std::vector<std::pair<double, std::pair<node, node> > > scores;
	G.forEdges([&](node u, node v) {
		double modScore = modScoring.edgeScore(u, v);
		std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
		scores.push_back(std::make_pair(-modScore, edge)); // negative due to minPQ
	});

	Aux::PriorityQueue<double, std::pair<node, node> > pq(scores);

	DEBUG("Initial PQ built before CNM loop");

// #define SLOW_PQ

	while (pq.size() > 0) {
		// determine best edge
		std::pair<double, std::pair<node, node> > pqMin = pq.extractMin();
		node best_u = none;
		node best_v = none;

#ifdef RECALC_ALL
		scores.clear();
		G.forEdges([&](node u, node v) {
			if (u != v) {
				double modScore = modScoring.edgeScore(u, v);
				std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
				scores.push_back(std::make_pair(-modScore, edge)); // negative due to minPQ
			}
		});
		Aux::PriorityQueue<double, std::pair<node, node> > pq2(scores);

		if (pq.size() == 0) {
			break;
		}

		std::pair<double, std::pair<node, node> > pqMin2 = pq2.extractMin();

		double bestDelta = std::numeric_limits<double>::lowest();

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

		if ((bestDelta != -pqMin.first) || (pqMin2.first != pqMin.first)) {
			pq.print();
			pq2.print();

			INFO("bestDelta: " , bestDelta , " at edge " , best_u , "," , best_v);
			best_u = (pqMin.second).first;
			best_v = (pqMin.second).second;
			INFO("pqMin.first: " , pqMin.first , " at edge " , best_u , "," , best_v , ", pqMin2.first: " , pqMin2.first);

			assert(pqMin.first == -bestDelta);
			assert(pqMin.first == pqMin2.first);
		}
#endif

		best_u = (pqMin.second).first;
		best_v = (pqMin.second).second;

		DEBUG("best_u: " , best_u , ", best_v: " , best_v);
		TRACE("rewiring of edges incident to best_v and removing them from PQ");

		// remove best_u edges from PQ (has to be done first)
		G.forWeightedEdgesOf(best_u, [&](node best_u, node neighbor, edgeweight ew) {
			std::pair<node, node> edge = (best_u <= neighbor) ? (std::make_pair(best_u, neighbor)) : (std::make_pair(neighbor, best_u));
			if (best_u != neighbor) {
				DEBUG("remove from PQ: " , -ew , " for " , best_u , " and " , neighbor);
				pq.remove(std::make_pair(ew, edge));
			}
		});

		// remove best_v edges from PQ (has to be done before rewiring)
		G.forWeightedEdgesOf(best_v, [&](node best_v, node neighbor, edgeweight ew) {
			std::pair<node, node> edge = (best_v <= neighbor) ? (std::make_pair(best_v, neighbor)) : (std::make_pair(neighbor, best_v));
			if (best_v != neighbor) {
				DEBUG("remove from PQ: " , -ew , " for " , best_v , " and " , neighbor);
				pq.remove(std::make_pair(ew, edge));
			}
		});


		// rewiring of edges incident to best_v
		G.forWeightedEdgesOf(best_v, [&](node best_v, node neighbor, edgeweight ew) {
			node endPoint = (best_v == neighbor) ? (best_u) : (neighbor);
			if (G.hasEdge(best_u, endPoint)) {
				// increase weight of existing edge
				G.increaseWeight(best_u, endPoint, ew);
				DEBUG("increase weight of edge " , best_u , "," , endPoint , " by " , ew);
			}
			else {
				// insert new edge to best_u
				G.addEdge(best_u, endPoint, ew);
				DEBUG("insert edge " , best_u , "," , endPoint);
			}
		});

		TRACE("remove best_v from graph");

		// remove best_v from graph
		G.forEdgesOf(best_v, [&](node best_v, node neighbor) {
			G.removeEdge(best_v, neighbor);
		});
		G.removeNode(best_v);

		// merge clusters best_u and best_v
		clustering.mergeClusters(clustering.clusterOf(best_u), clustering.clusterOf(best_v));

		// update delta mod scores for incident edges and insert them into the PQ, needs to
		// done after all edges have been rewired
#ifdef SLOW_PQ
		pq.clear();
		G.forEdges([&](node best_u, node neighbor) {
#else
		G.forEdgesOf(best_u, [&](node best_u, node neighbor) {
#endif
			if (best_u != neighbor) {
				std::pair<node, node> edge = (best_u <= neighbor) ? (std::make_pair(best_u, neighbor)) : (std::make_pair(neighbor, best_u));
				double deltaMod = modScoring.edgeScore(best_u, neighbor);
				DEBUG("insert into PQ: " , -deltaMod , " for " , best_u , " and " , neighbor);
				std::pair<double, std::pair<node, node> > elem = std::make_pair(-deltaMod, edge);
				pq.insert(elem);
			}
		});

		// compute new modularity
		double newModularity = modularityInspector.getQuality(clustering, graph);

		// record best solutions
		if (newModularity > bestModularity) {
			TRACE("new mod value: " , newModularity);
			bestModularity = newModularity;
			bestClustering = clustering;
		}
	}

	return bestClustering;
}

} // namespace

