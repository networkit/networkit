/*
 * LPDegreeOrdered.h
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#ifndef LPDEGREEORDERED_H_
#define LPDEGREEORDERED_H_

#include "CommunityDetectionAlgorithm.h"

namespace NetworKit {

typedef index label; // a label is the same as a cluster id

/**
 * @ingroup community
 * Label propagation-based community detection algorithm which
 * processes nodes in increasing order of node degree.
 */
class LPDegreeOrdered: public NetworKit::CommunityDetectionAlgorithm {
private:
	count nIterations = 0;	//!< number of iterations in last run


public:
	/**
	 * Constructor to the degree ordered label propagation community detection algorithm.
	 *
	 * @param[in]	G	input graph
	 */
	LPDegreeOrdered(const Graph& G);

	/**
	 * Detect communities.
	 */
	virtual void run() override;

	/**
	* Get number of iterations in last run.
	*
	* @return Number of iterations.
	*/
	virtual count numberOfIterations();

	std::string toString() const override;

};

} /* namespace NetworKit */
#endif /* LPDEGREEORDERED_H_ */
