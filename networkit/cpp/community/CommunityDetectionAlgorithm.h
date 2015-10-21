/*
 * CommunityDetectionAlgorithm.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef COMMUNITYDETECTIONALGORITHM_H_
#define COMMUNITYDETECTIONALGORITHM_H_

#include "../structures/Partition.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for community detection/graph clustering algorithms.
 */
class CommunityDetectionAlgorithm : public Algorithm {
public:
	/**
	 * A community detection algorithm operates on a graph, so the constructor expects a graph.
	 *
	 * @param[in]	G	input graph
	 */
	CommunityDetectionAlgorithm(const Graph& G);
	
	/**
	 * A community detection algorithm operates on a graph, so the constructor expects a graph.
	 *
	 * @param[in]	G	input graph
	 * @param[in]	baseClustering optional; the algorithm will start from the given clustering.
	 */
	CommunityDetectionAlgorithm(const Graph& G, const Partition baseClustering);

	/** Default destructor */
	virtual ~CommunityDetectionAlgorithm() = default;

	/**
	 * Apply algorithm to graph
	 */
	virtual void run() = 0;

	/**
	 * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
	 * @return partition of the node set
	 */
	virtual Partition getPartition();

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;

protected:
	const Graph& G;
	Partition result;
};

} /* namespace NetworKit */
#endif // COMMUNITYDETECTIONALGORITHM_H_
