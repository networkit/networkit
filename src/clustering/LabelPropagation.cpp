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

Clustering& LabelPropagation::run(Graph& G) {

	// TODO: remove after debugging
	using Aux::operator<<; // map printer

	typedef cluster label;	// a label is the same as a cluster id

	int64_t n = G.numberOfNodes();

	// create the clustering to be returned
	// set unique label for each node
	Clustering* labels = new Clustering(n);
	G.forallNodes([&](node v){
		labels->toSingleton(v);
	});

	int64_t majorityLabelCount = 0;	// number of nodes which already have the majority label
	int64_t nIterations = 0; 	// number of iterations

	// propagate labels
	while (majorityLabelCount != n) {
		majorityLabelCount = 0;
		nIterations += 1;


		// TODO: in random order
		G.forallNodes([&](node v){

			int64_t degV = G.degree(v);
			// ignore isolated nodes TODO: correct?
			if (degV > 0) {

				std::map<label, int64_t> labelCounts; // neighborLabelCounts maps label -> frequency in the neighbors

				std::cout << "new labelCounts: " << labelCounts << std::endl;

				// count the labels in the neighborhood of v and select the most frequent one
				G.forallNeighborsOf(v, [&](node w) {
					label cw = labels->clusterOf(w);
					if (labelCounts.find(cw) == labelCounts.end()) {
						labelCounts[cw] = 0;
					}
					labelCounts[cw] += 1;
				});

				std::cout << "labelCounts: " << labelCounts << std::endl;


				// get most frequent label
				label mostFrequent = 0; // TODO: check if 0 occurs in final clustering
				int64_t max = 0;
				for (auto it = labelCounts.begin(); it != labelCounts.end(); it++) {
					if (it->second > max) {
						max = it->second;
						mostFrequent = it->first;
					}
				}

				labels[v] = mostFrequent;
				// stop if v has label of at least half of its neighbors
				label dominantLabel = 0; // = None
				// try to find dominant label
				DEBUG("labelCounts size: " << labelCounts.size());
				for (auto it2 = labelCounts.begin(); it2 != labelCounts.end(); it2++) {
						DEBUG("labelCounts entry: " << it2->first << ":" << it2->second);
						if (it2->second > (degV / 2.0)) {
						dominantLabel = it2->first;
					}
				}

				if (dominantLabel != 0) {
					if ((*labels)[v] == dominantLabel) {
						majorityLabelCount += 1;
					}
				} // if no label dominant, do nothing
			}

		});

	}

	assert (labels != NULL);
	return *labels;

}

} /* namespace EnsembleClustering */
