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
			// count the labels in the neighborhood of v and select the most frequent one
			int64_t degV = G.degree(v);
			IndexMap<label, int64_t> neighborLabelCounts(n); // neighborLabelCounts[v] maps label -> frequency in the neighbors of v

			// TODO: forall neighbors of v
			//		count labels

			label mostFrequent;


			labels[v] = mostFrequent;
			// stop if v has label of at least half of its neighbors
			label dominantLabel = 0; // = None

			if (dominantLabel != 0) {
				if ((*labels)[v] == dominantLabel) {
					majorityLabelCount += 1;
				}
			} // if no label dominant, do nothing

		});

	}

	assert (labels != NULL);
	return *labels;

}

} /* namespace EnsembleClustering */
