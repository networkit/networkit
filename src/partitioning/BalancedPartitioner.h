/*
 * BalancedPartitioner.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef BALANCEDPARTITIONER_H_
#define BALANCEDPARTITIONER_H_

#include "../clustering/Clustering.h"
#include "../clustering/ClusteringGenerator.h"
#include "../clustering/EdgeCut.h"
#include "../io/ClusteringWriter.h"
#include "../clustering/Clustering.h"
#include "../auxiliary/RandomProbability.h"
#include "../coarsening/MatchingContracter.h"
#include "../matching/ParallelMatcher.h"
#include "../coarsening/ClusteringProjector.h"

namespace NetworKit {

class BalancedPartitioner {
public:
	BalancedPartitioner();
	virtual ~BalancedPartitioner();

	virtual Clustering run(Graph& G, count numBlocks) = 0;
	virtual Clustering& rerun(Graph& G, count numBlocks, Clustering& partition) = 0;

	virtual Clustering multilevelRun(Graph& graph, count numParts);
	virtual Clustering& multilevelRerun(Graph& graph, count numParts, Clustering& partition);

	virtual Clustering& postsmooth(Graph& graph, count numBlocks, Clustering& partition) = 0;
};

} /* namespace NetworKit */
#endif /* BALANCEDPARTITIONER_H_ */
