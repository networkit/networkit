/*
 * SelectiveCommunityDetector.h
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef SELECTIVECOMMUNITYDETECTOR_H_
#define SELECTIVECOMMUNITYDETECTOR_H_

#include <unordered_set>

#include "../auxiliary/Timer.h"
#include "../graph/Graph.h"

namespace NetworKit {

class SelectiveCommunityDetector {

public:

	SelectiveCommunityDetector(const Graph& G);

    /**
     * Detect communities for given seed nodes.
     * @return a mapping from seed node to community (as a set of nodes)
     */
	virtual std::map<node, std::set<node> >  run(std::set<unsigned int>& seeds) = 0;

	// FIXME: resolve Cython issue that does not allow a uint64_t as content type of a container as input


protected:

	const Graph& G;	//!< the input graph
};

} /* namespace NetworKit */
#endif /* SELECTIVECOMMUNITYDETECTOR_H_ */
