/*
 * Clusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef COMMUNITYDETECTIONALGORITHM_H_
#define COMMUNITYDETECTIONALGORITHM_H_

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for community detection/graph clustering algorithms.
 */
class CommunityDetectionAlgorithm {
public:
	/**
	 * A community detection algorithm operates on a graph, so the constructor expects a graph.
	 *
	 * @param[in]	G	input graph
	 */
	CommunityDetectionAlgorithm(const Graph& G);

	/** Default destructor */
	virtual ~CommunityDetectionAlgorithm() = default;

	/**
	 * Apply algorithm to graph
	 * @return partition of the node set
	 */
	virtual Partition run() = 0;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;

protected:
	const Graph& G;
};

} /* namespace NetworKit */
#endif /* CLUSTERER_H_ */
