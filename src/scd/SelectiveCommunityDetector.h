/*
 * SelectiveCommunityDetector.h
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef SELECTIVECOMMUNITYDETECTOR_H_
#define SELECTIVECOMMUNITYDETECTOR_H_

#include <unordered_set>

#include "../graph/Graph.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

class SelectiveCommunityDetector {

public:

	SelectiveCommunityDetector();

	virtual ~SelectiveCommunityDetector();

	virtual std::unordered_map<node, std::unordered_set<node> > run(Graph& G, std::unordered_set<node> seeds) = 0;
};

} /* namespace NetworKit */
#endif /* SELECTIVECOMMUNITYDETECTOR_H_ */
