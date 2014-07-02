/*
 * PLP.cpp
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "PLP.h"

#include <omp.h>
#include "../Globals.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/ProgressMeter.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/Random.h"
#include "../graph/NodeMap.h"

namespace NetworKit {

PLP::PLP(count theta) : updateThreshold(theta) {

	this->VERSION = "1.0";
}


Partition& PLP::runFromGiven(Graph& G, Partition& labels) {
	typedef index label; // a label is the same as a cluster id

	count n = G.numberOfNodes();
	index z = G.upperNodeIdBound();
	// update threshold heuristic
	if (updateThreshold == none) {
		updateThreshold = (count) (n / 1e5);
	}

	count nUpdated; // number of nodes which have been updated in last iteration
	nUpdated = n; // all nodes have new labels -> first loop iteration runs

	nIterations = 0; // number of iterations

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

	std::vector<bool> activeNodes(z); // record if node must be processed
	activeNodes.assign(z, true);

	Aux::Timer runtime;

	// propagate labels
	while (nUpdated > this->updateThreshold) { // as long as a label has changed...
		runtime.start();
		nIterations += 1;
		INFO("[BEGIN] LabelPropagation: iteration #" , nIterations);

		// reset updated
		nUpdated = 0;

		G.balancedParallelForNodes([&](node v){
			if ((activeNodes[v]) && (G.degree(v) > 0)) {

				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors

				// weigh the labels in the neighborhood of v
				G.forWeightedNeighborsOf(v, [&](node w, edgeweight weight) {
					label lw = labels.subsetOf(w);
					labelWeights[lw] += weight; // add weight of edge {v, w}
				});

				// get heaviest label
				label heaviest = std::max_element(labelWeights.begin(),
								labelWeights.end(),
								[](const std::pair<label, edgeweight>& p1, const std::pair<label, edgeweight>& p2) {
									return p1.second < p2.second;})->first;

				if (labels.subsetOf(v) != heaviest) { // UPDATE
					labels.moveToSubset(heaviest,v); //labels[v] = heaviest;
					nUpdated += 1; // TODO: atomic update?
					G.forNeighborsOf(v, [&](node u) {
						activeNodes[u] = true;
					});
				} else {
					activeNodes[v] = false;
				}

			} else {
				// node is isolated
			}
		});

		// for each while loop iteration...

		runtime.stop();
		DEBUG("[DONE] LabelPropagation: iteration #" , nIterations , " - updated " , nUpdated , " labels, time spent: " , runtime.elapsedTag());


	} // end while

	return labels;
}


Partition PLP::run(Graph& G) {
	// set unique label for each node
	index z = G.upperNodeIdBound();
	Partition labels(z);
	labels.allToSingletons();
	// TODO: make (call to) allToSingletons faster
//	G.parallelForNodes([&](node v) {
//		labels[v] = v;
//	});
//	labels.setUpperBound(z);

	return runFromGiven(G, labels);
}

std::string PLP::toString() const {
	std::stringstream strm;
	strm << "PLP(updateThreshold=" << this->updateThreshold << ")";
	return strm.str();
}


void PLP::setUpdateThreshold(count th) {
	this->updateThreshold = th;
}


count PLP::numberOfIterations() {
	return this->nIterations;
}

} /* namespace NetworKit */

