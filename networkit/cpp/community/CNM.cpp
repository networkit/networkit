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
#include <sstream>
#include "../auxiliary/Log.h"

namespace NetworKit{
CNM::CNM(const Graph& G) : CommunityDetectionAlgorithm(G) {};

node CNM::mergeEdge(Graph &G, node u, node v, bool discardSelfLoop){

	DEBUG("merge edge with nodes ", u," and ", v);

	if (u != v) {
		node newNode = G.addNode();

		// self-loop if necessary
		if (! discardSelfLoop) {
			TRACE("selfLoopWeight");
			edgeweight selfLoopWeight = G.weight(u, u) + G.weight(v, v) + G.weight(u, v);
			G.addEdge(newNode, newNode, selfLoopWeight);
			TRACE("end selfLoopWeight");
		}

		// rewire edges from u to newNode
		G.forEdgesOf(u, [&](node u, node neighbor, edgeweight w) {
			if (neighbor != u) {
				TRACE("neighbor of ",u,": ",neighbor);
				G.addEdge(neighbor, newNode, G.weight(u, neighbor)); // TODO: make faster
				TRACE("end neighbor of u");
			}
		});

		// rewire edges from v to newNode
		G.forEdgesOf(v, [&](node v, node neighbor, edgeweight w) {
			if (neighbor != v) {
				TRACE("neighbor of ",v,": ",neighbor);
				G.addEdge(neighbor, newNode, G.weight(v, neighbor));  // TODO: make faster
				TRACE("end neighbor of v");
			}
		});

		// delete edges of nodes to delete
		G.forEdgesOf(u, [&](node u, node neighbor) {
			G.removeEdge(u, neighbor);
		});
		G.forEdgesOf(v, [&](node v, node neighbor) {
			G.removeEdge(v, neighbor);
		});

		// delete nodes
		G.removeNode(u);
		G.removeNode(v);

		return newNode;
	}

	// no new node created
	return none;
}

void CNM::run() {
	// copy graph because we make changes due to merges
	Graph Gcopy(G.numberOfNodes(), true); // make weighted copy
	G.forEdges([&](node u, node v, edgeweight w){
		Gcopy.addEdge(u, v, w);
	});
	count n = Gcopy.numberOfNodes();

	auto createEdge([&](node u, node v) {
		std::pair<node, node> edge = (u <= v) ? (std::make_pair(u, v)) : (std::make_pair(v, u));
		return edge;
	});


	// start with singleton clustering
	Partition clustering(n);
	clustering.allToSingletons();

	// record current modularity
	double gTotalEdgeWeight = Gcopy.totalEdgeWeight();
	Modularity modularityInspector;
	modularityInspector.setTotalEdgeWeight(gTotalEdgeWeight);
	double bestModularity = modularityInspector.getQuality(clustering, Gcopy);
	Partition bestClustering = clustering;

	// insert edges into priority queue, key: delta mod score
	ModularityScoring<double> modScoring(Gcopy);
	std::vector<double> scores;
	std::map<std::pair<node, node>, index> edgeToIndex;
	std::map<index, std::pair<node, node> > indexToEdge;
	index idx = 0;
	Gcopy.forEdges([&](node u, node v) {
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
//		TRACE("best_u: " , best_u , ", best_v: " , best_v, ", ew: ", Gcopy.weight(best_u, best_v), ", deltaMod: ", -pqOpt.first);

		// merge clusters best_u and best_v
		index newCluster = clustering.mergeSubsets(clustering.subsetOf(best_u), clustering.subsetOf(best_v));

		// adapt PQ accordingly
		Gcopy.forEdgesOf(best_u, [&](node best_u, node neighbor) {
			pq.remove(edgeToIndex[createEdge(best_u, neighbor)]);
		});
		Gcopy.forEdgesOf(best_v, [&](node best_v, node neighbor) {
			pq.remove(edgeToIndex[createEdge(best_v, neighbor)]);
		});

		// merge incident nodes in Gcopy
		node newNode = CNM::mergeEdge(Gcopy, best_u, best_v, false);
		assert(newNode != none);

		// adapt clustering accordingly
		clustering.extend();
		clustering[newNode] = newCluster;

		// update delta mod scores for incident edges and insert them into the PQ,
		// needs to done after all edges have been rewired
		Gcopy.forEdgesOf(newNode, [&](node newNode, node neighbor) {
			if (newNode != neighbor) {
				// determine edge score
//				assert(Gcopy.totalEdgeWeight() == gTotalEdgeWeight);
				ModularityScoring<double> modScoring(Gcopy, gTotalEdgeWeight);
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
//		assert(Gcopy.totalEdgeWeight() == gTotalEdgeWeight);
		modularityInspector.setTotalEdgeWeight(gTotalEdgeWeight);
		double newModularity = modularityInspector.getQuality(clustering, G);
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
	result = std::move(bestClustering);
	hasRun = true;
}

} // namespace
