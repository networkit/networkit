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

	int64_t nUpdated;	// number of nodes which have been updated in last iteration
	nUpdated = n; // all nodes have new labels -> first loop iteration runs

	int64_t nIterations = 0; 	// number of iterations
	int64_t nDominated = 0;	// number of nodes which are already dominated



	/**
	 * == Dealing with isolated nodes ==
	 *
	 * The pseudocode published does not deal with isolated nodes (and therefore does not terminate if they are present).
	 * Isolated nodes stay singletons. They can be ignored in the while loop, but the loop condition must
	 * compare to the number of non-isolated nodes instead of n.
	 *
	 * == Termination criterion ==
	 *
	 * The published termination criterion is: All nodes have got the label of the majority of their neighbors.
	 * In general this does not work. It was changed to: No label was changed in last iteration.
	 */


	// PERFORMANCE: precompute and store incident edge weight for all nodes
	INFO("[BEGIN] Label Propagation: precomputing incident weight");
	Aux::ProgressMeter pm(n, 1000);
	NodeMap<double> incidentWeight(n, 0.0);
	G.forallNodes([&](node v) {
		incidentWeight[v] = G.incidentWeight(v);
		pm.signal(v);
	}, "parallel");
	pm.end();


	// propagate labels
	while (nUpdated > 0) { // as long as a label has changed...
		nIterations += 1;
		INFO("[BEGIN] LabelPropagation: iteration #" << nIterations);

		// reset updated
		nUpdated = 0;


		NodeMap<int> dominated(n, 0); // map node-> dominated?

		std::vector<node> shuffledNodes;
		shuffledNodes.resize(n); 	// hold n nodes
		G.forallNodes([&](node v){
			shuffledNodes[v - 1] = v;	// store all nodes in vector
		});
		std::shuffle(shuffledNodes.begin(), shuffledNodes.end(), randgen);

		Aux::ProgressMeter pm(shuffledNodes.size(), 1000);
		int64_t vc = 0;  // node counter

		for (node v : shuffledNodes) {
			// PROGRESS
			vc += 1;
			pm.signal(vc);


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
				label heaviest = 0;
				double maxWeight = 0.0;
				for (auto it = labelWeights.begin(); it != labelWeights.end(); it++) {
					if (it->second > maxWeight) {
						maxWeight = it->second;
						heaviest = it->first;
					}
				}

				if (labels[v] != heaviest) { // UPDATE
					// DEBUG
					TRACE("updating label of " << v << " from " << labels[v] << " to " << heaviest);
					// DEBUG
					labels[v] = heaviest;
					nUpdated += 1;
				} else {
					TRACE("label of " << v << " stays " << labels[v]);
				}

			} else {
				// node is isolated
				TRACE("ignoring isolated node: " << v);
			}

		} // end for shuffled nodes

		// PROGRESS
		pm.end();
		// for each while loop iteration...

	} // end while

	return labels;

}

} /* namespace EnsembleClustering */
