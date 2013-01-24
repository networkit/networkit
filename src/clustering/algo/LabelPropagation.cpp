/*
 * LabelPropagation.cpp
 *
 *  Created on: 07.12.2012
 *      Author: cls
 */

#include "LabelPropagation.h"

namespace EnsembleClustering {

LabelPropagation::LabelPropagation() {
	// TODO Auto-generated constructor stub

}

LabelPropagation::~LabelPropagation() {
	// TODO Auto-generated destructor stub
}

Clustering LabelPropagation::run(Graph& G) {


	typedef cluster label;	// a label is the same as a cluster id

	// init random for std::shuffle
	std::default_random_engine rd;
	std::mt19937 randgen(rd());

	int64_t n = G.numberOfNodes();

	// create the clustering to be returned
	// set unique label for each node
	Clustering labels(n);
	labels.allToSingletons();

	int64_t majorityLabelCount = 0;	// number of nodes which already have the majority label
	int64_t nIterations = 0; 	// number of iterations

	/**
	 * == Dealing with isolated nodes ==
	 *
	 * The pseudocode published does not deal with isolated nodes (and therefore does not terminate if they are present).
	 * Isolated nodes stay singletons. They can be ignored in the while loop, but the loop condition must
	 * compare to the number of non-isolated nodes instead of n.
	 *
	 */

	// count connected nodes for loop condition
	int64_t nConnected = 0;
	G.forallNodes([&](node v) {
		if (G.degree(v) > 0) {
			#pragma omp atomic update
			nConnected += 1;
		}
	}, "parallel");

	// PERFORMANCE: precompute and store incident edge weight for all nodes
	NodeMap<double> incidentWeight(n);
	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
	}, "parallel");


	// propagate labels
	while (majorityLabelCount != nConnected) {
		nIterations += 1;
		DEBUG("***** LabelPropagation: iteration #" << nIterations << "*****");
		// DEBUG
		TRACE("number of nodes which already have the majority label: " << majorityLabelCount << " of " << G.numberOfNodes());
		// DEBUG

		// reset majority label count
		majorityLabelCount = 0;
		NodeMap<int> dominated(nConnected);

		// DEBUG
		if (nIterations >= 42) {
			ERROR("LabelPropagation reached " << nIterations << " iterations. It usually terminates after less than 5 iterations. Something has gone terribly wrong.");
			ERROR("majorityLabelCount=" << majorityLabelCount);
			GraphIO graphio;
			graphio.writeAdjacencyList(G, "sandbox/LabelPropagationFAIL.adjlist");
			throw std::runtime_error("aborting LabelPropagation to avoid infinite loop");
		}
		// DEBUG

		std::vector<node> shuffledNodes;
		shuffledNodes.resize(n); 	// hold n nodes
		G.forallNodes([&](node v){
			shuffledNodes[v - 1] = v;	// store all nodes in vector
		}, "parallel");
		std::shuffle(shuffledNodes.begin(), shuffledNodes.end(), randgen);
		// DEBUG
		TRACE("shuffledNodes: " << Aux::vectorToString(shuffledNodes));
		// DEBUG

		for (node v : shuffledNodes) {
			// ignore isolated nodes TODO: correct?
			if (G.degree(v) > 0) {

				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors


				// weigh the labels in the neighborhood of v
				G.forallNeighborsOf(v, [&](node w) {
					label lw = labels[w];
					if (labelWeights.find(lw) == labelWeights.end()) {
						labelWeights[lw] = 0.0; // init map entry if not yet in map
					}
					labelWeights[lw] += G.weight(v, w);	// add weight of edge {v, w}
				});

				// consider also self-loop (i.e. v's own weight)
				label lv = labels[v];
				if (labelWeights.find(lv) == labelWeights.end()) {
					labelWeights[lv] = 0.0;	// init map entry if not yet in map
				}
				labelWeights[lv] += G.weight(v);


				// get most frequent label
				label heaviest = 0; // TODO: check if 0 occurs in final clustering
				double maxWeight = 0.0;
				for (auto it = labelWeights.begin(); it != labelWeights.end(); it++) {
					if (it->second > maxWeight) {
						maxWeight = it->second;
						heaviest = it->first;
					}
				}

				if (labels[v] != heaviest) {
					TRACE("updating label of " << v << " from " << labels[v] << " to " << heaviest);
					labels.moveToCluster(heaviest, v);
				} else {
					TRACE("label of " << v << " stays " << labels[v]);
				}

				// stop if v has label of at least half of its neighbors
				label dominantLabel = 0; // = None
				// try to find dominant label
				for (auto it2 = labelWeights.begin(); it2 != labelWeights.end(); it2++) {
					// label is dominant if it weighs more than half of the incident weight (including self-loop)
						if (it2->second > ((incidentWeight[v] + G.weight(v)) / 2.0)) {
						dominantLabel = it2->first;
					}
				}

				if (dominantLabel != 0) {
					if (labels[v] == dominantLabel) {
						majorityLabelCount += 1;
						dominated[v] = 1;
					}
				} else {
					TRACE("no dominant label found for node: " << v);
				}// if no label dominant, do nothing

			} else {
				// node is isolated
				TRACE("ignoring isolated node: " << v);
			}
		} // end for shuffled nodes

		DEBUG("dominated: " << dominated.toString());
	} // end while

	return labels;

}

} /* namespace EnsembleClustering */
