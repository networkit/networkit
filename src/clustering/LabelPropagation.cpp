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

	typedef cluster label;	//!< a label is the same as a cluster id

	int64_t n = G.numberOfNodes();

	std::unordered_map<label, int64_t> neighborLabelCounts[n]; //!< neighborLabelCounts[v] maps label -> frequency in the neighbors of v

	Clustering* labels = new Clustering(n);

	node v;

	// TODO: for all nodes
	labels->toSingleton(v);

	int64_t majorityLabelCount = 0;	//!< number of nodes which already have the majority label
	int64_t nIterations = 0; 	//!< number of iterations

	// TODO: for all nodes
	int64_t degV = G.getDegree(v);
	neighborLabelCounts[v].clear();
		// TODO: for all neighbors w
	node w;

	if (neighborLabelCounts[v].count(labels[w]) == 0) {
		neighborLabelCounts[v][labels[w]] = 1;
	} else {
		neighborLabelCounts[v][labels[w]]++;
	}

	assert (labels != NULL);
	return *labels;

}

} /* namespace EnsembleClustering */
