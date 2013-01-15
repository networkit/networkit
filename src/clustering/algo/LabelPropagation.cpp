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

	// TODO: remove after debugging
	using Aux::operator<<; // map printer

	typedef cluster label;	// a label is the same as a cluster id

	int64_t n = G.numberOfNodes();

	// create the clustering to be returned
	// set unique label for each node
	Clustering labels(n);
	// TODO: use labels.allToSingleton
	G.forallNodes([&](node v){
		labels.toSingleton(v);
	});

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
			nConnected += 1;
		}
	});


	// propagate labels
	while (majorityLabelCount != nConnected) {
		nIterations += 1;
		DEBUG("***** LabelPropagation: iteration #" << nIterations << "*****");
		// DEBUG
		TRACE("number of nodes which already have the majority label: " << majorityLabelCount << " of " << G.numberOfNodes());
		// DEBUG

		// reset majority label count
		majorityLabelCount = 0;

		// DEBUG
		if (nIterations >= 42) {
			ERROR("LabelPropagation reached " << nIterations << " iterations. It usually terminates after less than 5 iterations. Something has gone terribly wrong.");
			throw std::runtime_error("aborting LabelPropagation to avoid infinite loop");
		}
		// DEBUG

		std::vector<node> shuffledNodes;
		G.forallNodes([&](node v){
			shuffledNodes.push_back(v);
		});
		std::random_shuffle(shuffledNodes.begin(), shuffledNodes.end());
		TRACE("shuffledNodes: " << Aux::vectorToString(shuffledNodes));

		for (node v : shuffledNodes) {
			// ignore isolated nodes TODO: correct?
			if (G.degree(v) > 0) {

				std::map<label, int64_t> labelCounts; // neighborLabelCounts maps label -> frequency in the neighbors


				// count the labels in the neighborhood of v and select the most frequent one
				G.forallNeighborsOf(v, [&](node w) {
					label cw = labels.clusterOf(w);
					if (labelCounts.find(cw) == labelCounts.end()) {
						labelCounts[cw] = 0;
					}
					labelCounts[cw] += 1;
				});

				// get most frequent label
				label mostFrequent = 0; // TODO: check if 0 occurs in final clustering
				int64_t max = 0;
				for (auto it = labelCounts.begin(); it != labelCounts.end(); it++) {
					if (it->second > max) {
						max = it->second;
						mostFrequent = it->first;
					}
				}

				TRACE("updating label of " << v << " from " << labels.clusterOf(v) << " to " << mostFrequent);
				labels.moveToCluster(mostFrequent, v);

				// stop if v has label of at least half of its neighbors
				label dominantLabel = 0; // = None
				// try to find dominant label
				for (auto it2 = labelCounts.begin(); it2 != labelCounts.end(); it2++) {
						if (it2->second > (G.degree(v) / 2.0)) {
						dominantLabel = it2->first;
					}
				}

				if (dominantLabel != 0) {
					if (labels[v] == dominantLabel) {
						majorityLabelCount += 1;
					}
				} else {
					TRACE("no dominant label found for node: " << v);
				}// if no label dominant, do nothing
			} else {
				// node is isolated
				TRACE("ignoring isolated node: " << v);
			}
		} // end for shuffled nodes

	}

	return labels;

}

} /* namespace EnsembleClustering */
