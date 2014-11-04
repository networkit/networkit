/*
 * CNM.h
 *
 *  Created on: Jun 10, 2013
 *      Author: michi
 */

#ifndef CNM_H_
#define CNM_H_

#include "CommunityDetectionAlgorithm.h"
#include <utility>
#include <stdexcept>
#include <utility>
#include "../auxiliary/Log.h"

namespace NetworKit {

/**
 * @ingroup community
 * Clustering algorithm due to Clauset, Newman and Moore.
 * Probably not the fastest possible implementation, but it already uses a priority queue
 * and local updates.
 */
class CNM : public NetworKit::CommunityDetectionAlgorithm {
public:
	/**
	 * Constructor to the community detection algorithm after Clauset, Newman and Moore.
	 *
	 * @param[in]	G	input graph
	 */
	CNM(const Graph& G);

	/**
	 * Detect communities in the given graph @a graph.
	 *
	 * @return A partition containing the found communities.
	 */
	Partition run() override;

	/**
	 * Get string representation.
	 *
	 * @return A string representation of this algorithm.
	 */
	std::string toString() const override {
		return "CNM";
	}

protected:
	static node mergeEdge(Graph &G, node u, node v ,bool discardSelfLoop);

};

}

#endif /* CNM_H_ */
