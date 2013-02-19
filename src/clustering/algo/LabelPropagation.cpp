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


static int countOne(const bool& a) {
	return (int) a;
}

Clustering LabelPropagation::run(Graph& G) {
	typedef cluster label; // a label is the same as a cluster id

	// get global variables
	const bool printProgress = PRINT_PROGRESS;
	const bool randOrder = RAND_ORDER;

	// init random for std::shuffle
	std::default_random_engine rd;
	std::mt19937 randgen(rd());

	// open file for csv output
	std::stringstream filePath;
	filePath << "output/LPCount-" << G.getName() << ".csv";
	std::ofstream lpCount(filePath.str());
	lpCount << "nActive;nUpdated" << std::endl; // header

	count n = G.numberOfNodes();

	// create the clustering to be returned
	// set unique label for each node
	Clustering labels(n);
	labels.allToSingletons();

	count nUpdated; // number of nodes which have been updated in last iteration
	nUpdated = n; // all nodes have new labels -> first loop iteration runs

	count nIterations = 0; // number of iterations

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

	std::vector<node> shuffledNodes(n);
	G.parallelForNodes([&](node v) {
		shuffledNodes[v] = v; // store all nodes in vector
		});

	std::vector<bool> activeNodes(n); // record if node must be processed
	activeNodes.assign(n, true);

	Aux::Timer runtime;

	// propagate labels
	while (nUpdated > this->updateThreshold) { // as long as a label has changed...
		runtime.start();
		nIterations += 1;
		INFO("[BEGIN] LabelPropagation: iteration #" << nIterations);

		// reset updated
		nUpdated = 0;

		if (randOrder) {
			// new random order
#ifdef _OPENMP
			// not really random, but the next best thing in parallel w/o hurting performance
			count numChunks = omp_get_num_threads();
			count chunkSize = (count) n / numChunks;// discard remainder
#pragma omp parallel for
			for (index i = 0; i < numChunks; ++i) {
				index begin = i * chunkSize;
				index end = begin + chunkSize;
				std::shuffle(&shuffledNodes[begin], &shuffledNodes[end], randgen);
			}
#else
			std::shuffle(shuffledNodes.begin(), shuffledNodes.end(), randgen);
#endif
		}

		Aux::ProgressMeter pm(n, 10000);

		// TODO: delete for performance tests
		count nActive = std::count_if(activeNodes.begin(), activeNodes.end(), countOne);
		INFO("number of active nodes: " << nActive);


#pragma omp parallel for schedule(guided) shared(nUpdated)
		for (int64_t i = 0; i < n; ++i) {
			node v = shuffledNodes[i];

			// PROGRESS

			if (printProgress) {
				pm.signal(i);
			}

			if ((activeNodes[v]) && (G.degree(v) > 0)) {

				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors

				// weigh the labels in the neighborhood of v
				G.forWeightedNeighborsOf(v, [&](node w, edgeweight weight) {
					label lw = labels[w];
					labelWeights[lw] += weight; // add weight of edge {v, w}
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
		if (printProgress) {
			pm.end();
		}

		runtime.stop();
		INFO("[DONE] LabelPropagation: iteration #" << nIterations << " - updated " << nUpdated << " labels, time spent: " << runtime.elapsedTag());

		// record nActive and nUpdated in csv file
		lpCount << nActive << ";" << nUpdated << std::endl;

	} // end while

	return labels;

}

std::string LabelPropagation::toString() {
	std::stringstream strm;
	strm << "LabelPropagation(randOrder=" << RAND_ORDER << ",updateThreshold=" << this->updateThreshold << ")";
	return strm.str();

}

void LabelPropagation::setUpdateThreshold(count th) {
	this->updateThreshold = th;
}

} /* namespace EnsembleClustering */
