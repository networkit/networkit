/*
 * BalancedPartitioner.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef BALANCEDPARTITIONER_H_
#define BALANCEDPARTITIONER_H_

#include "../structures/Partition.h"
#include "../community/ClusteringGenerator.h"
#include "../community/EdgeCut.h"
#include "../io/ClusteringWriter.h"
#include "../coarsening/MatchingContracter.h"
#include "../matching/ParallelMatcher.h"
#include "../coarsening/ClusteringProjector.h"

namespace NetworKit {

class BalancedPartitioner {
protected:
	float balance;

public:
	BalancedPartitioner();
	virtual ~BalancedPartitioner();

	virtual Partition run(Graph& G, count numBlocks) = 0;
	virtual Partition& rerun(Graph& G, count numBlocks, Partition& partition) = 0;

	virtual Partition multilevelRun(Graph& graph, count numParts);
	virtual Partition& multilevelRerun(Graph& graph, count numParts, Partition& partition);

	virtual Partition& postsmooth(Graph& graph, count numBlocks, Partition& partition) = 0;

	void setBalance(float balanceFactor);
};

} /* namespace NetworKit */
#endif /* BALANCEDPARTITIONER_H_ */
