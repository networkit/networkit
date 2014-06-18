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
 * Abstract base class for community detection/graph clustering algorithms.
 */
class CommunityDetectionAlgorithm {
public:

	/** Default destructor */
	virtual ~CommunityDetectionAlgorithm();

	/**
	 * Apply algorithm to graph
	 * @return partition of the node set
	 */
	virtual Partition run(Graph& G) = 0;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif /* CLUSTERER_H_ */
