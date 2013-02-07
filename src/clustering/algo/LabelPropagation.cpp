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

	const bool printProgress = PRINT_PROGRESS;

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
	NodeMap<double> weightedDegree(n, 0.0);
	G.parallelForNodes([&](node v) {
		weightedDegree[v] = G.weightedDegree(v);
		if (printProgress) {
			pm.signal(v);
		}
	});
	if (printProgress) {
		pm.end();
	}



	// propagate labels
	while (nUpdated > 0) { // as long as a label has changed...
		nIterations += 1;
		INFO("[BEGIN] LabelPropagation: iteration #" << nIterations);

		// reset updated
		nUpdated = 0;

		std::vector<node> shuffledNodes;
		shuffledNodes.resize(n); 	// hold n nodes
		G.parallelForNodes([&](node v){
			shuffledNodes[v] = v;	// store all nodes in vector
		});
		std::shuffle(shuffledNodes.begin(), shuffledNodes.end(), randgen);

		Aux::ProgressMeter pm(n, 1000);

		#pragma omp parallel for
		for (int64_t i = 0; i < n; ++i) {
			node v = shuffledNodes[i];

			// PROGRESS

			if (printProgress) {
				pm.signal(i);
			}


			if (G.degree(v) > 0) {

				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors


				// weigh the labels in the neighborhood of v
				G.forNeighborWeightsOf(v, [&](node w, edgeweight weight) {
					label lw = labels[w];
					labelWeights[lw] += weight;	// add weight of edge {v, w}
				});


				// get heaviest label
				label heaviest = std::max_element(labelWeights.begin(), labelWeights.end(),
				      [](const std::pair<label, edgeweight>& p1, const std::pair<label, edgeweight>& p2) {
				        return p1.second < p2.second; })->first;

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

		// for each while loop iteration...

		// PROGRESS
		if (printProgress) {
			pm.end();
		}


		INFO("[DONE] LabelPropagation: iteration #" << nIterations << " - updated " << nUpdated << " labels");


	} // end while

	return labels;

}

} /* namespace EnsembleClustering */
