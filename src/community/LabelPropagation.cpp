/*
 * LabelPropagation.cpp
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "LabelPropagation.h"

#include <omp.h>
#include "../Globals.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/ProgressMeter.h"
#include "../auxiliary/Timer.h"
#include "../auxiliary/RandomInteger.h"
#include "../graph/NodeMap.h"

namespace NetworKit {

LabelPropagation::LabelPropagation(count theta) : updateThreshold(theta) {

	this->VERSION = "1.0";
}

LabelPropagation::~LabelPropagation() {
	// TODO Auto-generated destructor stub
}


Clustering LabelPropagation::run(Graph& G) {
	typedef cluster label; // a label is the same as a cluster id

	// get global variables
	const bool printProgress = PRINT_PROGRESS;
	const bool randOrder = RAND_ORDER;							// explicitly randomize node order for each iteration
	const count inactiveSeeds = INACTIVE_SEEDS;					// number of seed nodes which are set inactive for first iteration
	const bool normalizeVotes = NORMALIZE_VOTES;
	const bool scaleClusterStrength = (SCALE_STRENGTH != 0.0);
	std::vector<double> scale;

	// init random for std::shuffle
	std::default_random_engine rd;
	std::mt19937 randgen(rd());

	// option not used
	if (normalizeVotes) {
		WARN("normalized votes turned off for undirected graphs");
	}

	count n = G.numberOfNodes();

	// set unique label for each node
	Clustering labels(n);
	labels.allToSingletons();

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

	// PERFORMANCE: precompute and store incident edge weight for all nodes
	DEBUG("[BEGIN] Label Propagation: precomputing incident weight");
	NodeMap<double> weightedDegree(n, 0.0);
	G.parallelForNodes([&](node v) {
		weightedDegree[v] = G.weightedDegree(v);
	});



	// prepare nodes for randomized iteration order
	std::vector<node> nodes(n);
	G.parallelForNodes([&](node v) {
		nodes[v] = v; // store all nodes in vector
	});

	std::vector<bool> activeNodes(n); // record if node must be processed
	activeNodes.assign(n, true);


	// randomize outcome by deactivating a very small number of nodes
	// TODO: make this an object attribute
	if (inactiveSeeds > 0) {
		Aux::RandomInteger randInt;
		for (count i = 0; i < inactiveSeeds; i++) {
			node u = randInt.generate(0, (n-1));
			activeNodes[u] = false;
		}
	}





	Aux::Timer runtime;

	// propagate labels
	while (nUpdated > this->updateThreshold) { // as long as a label has changed...
		runtime.start();
		nIterations += 1;
		INFO("[BEGIN] LabelPropagation: iteration #" << nIterations);

		// reset updated
		nUpdated = 0;

		if (randOrder) { // randomize the order of node iteration
			// new random order
#ifdef _OPENMP
			// not really random, but the next best thing in parallel w/o hurting performance
			count numChunks = omp_get_num_threads();
			count chunkSize = (count) n / numChunks;// discard remainder
#pragma omp parallel for
			for (index i = 0; i < numChunks; ++i) {
				index begin = i * chunkSize;
				index end = begin + chunkSize;
				std::shuffle(&nodes[begin], &nodes[end], randgen);
			}
#else
			std::shuffle(nodes.begin(), nodes.end(), randgen);
#endif
		}

		if (scaleClusterStrength) {
			// TODO: documentation?
			std::vector<count> clusterSizes = labels.clusterSizes();
			scale.resize(clusterSizes.size());
			INFO("Scaling cluster strengths with exponent " << SCALE_STRENGTH);
			for (index i = 0; i < clusterSizes.size(); ++i) {
				scale[i] = pow((double) clusterSizes[i], SCALE_STRENGTH);
			}
		}

		//Aux::ProgressMeter pm(n, 10000);
		// Reason: performance decreases for very big graphs

		// removed for performance reasons
		// count nActive = std::count_if(activeNodes.begin(), activeNodes.end(), countOne);
		// INFO("number of active nodes: " << nActive);

// TODO: make this forNodes loop
#pragma omp parallel for schedule(guided) shared(nUpdated)
		for (index i = 0; i < n; ++i) {
			node v = nodes[i];

			// PROGRESS
			/*if (printProgress) {
				pm.signal(i);
			}*/

			if ((activeNodes[v]) && (G.degree(v) > 0)) {

				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors

				// weigh the labels in the neighborhood of v
				G.forWeightedNeighborsOf(v, [&](node w, edgeweight weight) {
					label lw = labels[w];
					if (scaleClusterStrength) {
						labelWeights[lw] += weight / scale[labels[v]]; // add weight of edge {v, w}
					} else {
						labelWeights[lw] += weight; // add weight of edge {v, w}
					}
				});

				// get heaviest label
				label heaviest = std::max_element(labelWeights.begin(),
								labelWeights.end(),
								[](const std::pair<label, edgeweight>& p1, const std::pair<label, edgeweight>& p2) {
									return p1.second < p2.second;})->first;

				if (labels[v] != heaviest) { // UPDATE
					labels[v] = heaviest;
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

		} // end for shuffled nodes

		// for each while loop iteration...

		// PROGRESS
		/*if (printProgress) {
			pm.end();
		}*/

		runtime.stop();
		DEBUG("[DONE] LabelPropagation: iteration #" << nIterations << " - updated " << nUpdated << " labels, time spent: " << runtime.elapsedTag());

	} // end while

	return labels;

}

std::string LabelPropagation::toString() const {
	std::stringstream strm;
	// FIXME: report current number of threads
#ifdef _OPENMP_
	strm << "LabelPropagation(" << "version=" << this->VERSION << ", randOrder=" << RAND_ORDER << ",updateThreshold=" << this->updateThreshold << ",OpenMP)";
#else
	strm << "LabelPropagation(" << "version=" << this->VERSION << ", randOrder=" << RAND_ORDER << ",updateThreshold=" << this->updateThreshold << ",threads=1)";
#endif
	return strm.str();
}

void LabelPropagation::setUpdateThreshold(count th) {
	this->updateThreshold = th;
}


count LabelPropagation::numberOfIterations() {
	return this->nIterations;
}

} /* namespace NetworKit */
