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
	virtual ~SelectiveCommunityDetector() = default;

    /**
     * Detect communities for given seed nodes.
     * @return a mapping from seed node to community (as a set of nodes)
     */
	virtual std::map<node, std::set<node> >  run(const std::set<node>& seeds) = 0;

protected:

	const Graph& G;	//!< the input graph
};

} /* namespace NetworKit */
#endif /* SELECTIVECOMMUNITYDETECTOR_H_ */
