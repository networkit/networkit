/*
 * BalancedLabelPropagation.h
 *
 *  Created on: Jun 19, 2013
 *      Author: Henning
 */

#ifndef BALANCEDLABELPROPAGATION_H_
#define BALANCEDLABELPROPAGATION_H_

#include "../community/ClusteringGenerator.h"
#include "../community/EdgeCut.h"
#include "../io/PartitionWriter.h"
#include "../structures/Partition.h"
#include "../coarsening/MatchingContracter.h"
#include "../matching/ParallelMatcher.h"
#include "../coarsening/ClusteringProjector.h"
#include "BalancedPartitioner.h"
#include <vector>
#include <algorithm>

namespace NetworKit {

/**
 * TODO: class documentation
 */
class BalancedLabelPropagation: public BalancedPartitioner {
protected:
	double exponent;

public:
	BalancedLabelPropagation(double exponent);
	virtual ~BalancedLabelPropagation();

	virtual Partition run(Graph& graph, count numBlocks);
	virtual Partition& rerun(Graph& graph, count numBlocks, Partition& partition);
//	virtual Partition multilevelRun(Graph& graph, count numParts);
//	virtual Partition& multilevelRerun(Graph& graph, count numParts, Partition& partition);

	virtual Partition& postsmooth(Graph& graph, count numBlocks, Partition& partition);

	void setExponent(double exponent) {
		this->exponent = exponent;
	}
};

} /* namespace NetworKit */
#endif /* BALANCEDLABELPROPAGATION_H_ */
