/*
 * KernighanLin.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#include "KernighanLin.h"

namespace NetworKit {

KernighanLin::KernighanLin() {
	// TODO Auto-generated constructor stub

}

KernighanLin::~KernighanLin() {
	// TODO Auto-generated destructor stub
}

Clustering KernighanLin::run(Graph& G, count numBlocks) {
	assert(numBlocks == 2); // FIXME: so far only bipartitioning
	ClusteringGenerator gen;
	Clustering partition = gen.makeContinuousBalancedClustering(G, numBlocks);
	partition = this->rerun(G, numBlocks, partition);
	return partition;
}

Clustering& KernighanLin::rerun(Graph& G, count numBlocks,
		Clustering& partition) {

	// compute initial gains

	// perform loop:

	// move each vertex once

	// use best result found


	// TODO: implement KL




	return partition;
}

Clustering& KernighanLin::postsmooth(Graph& graph, count numBlocks,
		Clustering& partition) {
	// TODO: possibly something smarter (if necessary)
	return partition;
}

} /* namespace NetworKit */
