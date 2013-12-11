/*
 * Clusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERER_H_
#define CLUSTERER_H_

#include "../clustering/Clustering.h"

namespace NetworKit {

/**
 * Abstract base class for community detection/graph clustering algorithms.
 */
class Clusterer {
public:

	Clusterer();

	virtual ~Clusterer();

	/**
	 * Apply algorithm to graph
	 * @return partition of the node set
	 */
	virtual Clustering run(Graph& G) = 0;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif /* CLUSTERER_H_ */
