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
#include "../base/Parameters.h"

namespace NetworKit {

class SelectiveCommunityDetector {

public:

	SelectiveCommunityDetector(const Graph& G);

	virtual ~SelectiveCommunityDetector();

	virtual void run(std::unordered_set<node> seeds) = 0;

	/** 
	 * @return a mapping from seed node to community (as a set of nodes)
	 */
	virtual std::unordered_map<node, std::unordered_set<node> > getResult() = 0;

	/** 
	 * @return time in milliseconds spent on processing each seed node
	 */
	virtual std::unordered_map<node, double> getTimings() = 0;

public:

	const Graph& G;	//!< the input graph
};

} /* namespace NetworKit */
#endif /* SELECTIVECOMMUNITYDETECTOR_H_ */
