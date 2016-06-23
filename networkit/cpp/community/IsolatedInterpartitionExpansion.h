#ifndef ISOLATEDINTERPARTITIONEXPANSION_H
#define ISOLATEDINTERPARTITIONEXPANSION_H

#include "LocalPartitionEvaluation.h"

namespace NetworKit {

/**
 * Isolated inter-partition expansion is a measure for how well a partition
 * (communtiy/cluster) is separated from the rest of the graph.
 *
 * The expansion of a partition is defined as the weight of the cut divided
 * by number of nodes in the partition or in the rest of the graph, whatever
 * is smaller. Small values thus indicate that the cut is small compared to
 * the size of the smaller of the separated parts. For the whole partitions
 * usually the maximum or the unweighted average is used. Note that expansion
 * values can be larger than 1.
 *
 * See also Experiments on Density-Constrained Graph Clustering,
 * Robert Görke, Andrea Kappes and Dorothea Wagner, JEA 2015:
 * http://dx.doi.org/10.1145/2638551
 */
class IsolatedInterpartitionExpansion : public LocalPartitionEvaluation {
public:
	using LocalPartitionEvaluation::LocalPartitionEvaluation;

	/**
	 * Execute the algorithm.
	 */
	virtual void run() override;

	/**
	 * @return true - smaller values are better than larger values.
	 */
	virtual bool isSmallBetter() const override { return true; };

	/**
	 * @return false - this algorithm is not parallel.
	 */
	virtual bool isParallel() const override { return false; };

	/**
	 * Get the name of the algorithm.
	 */
	virtual std::string toString() const override { return "Isolated inter-partition expansion"; };
};

}

#endif // ISOLATEDINTERPARTITIONEXPANSION_H
