/*
 * SelectiveSCAN.h
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#ifndef SELECTIVESCAN_H_
#define SELECTIVESCAN_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

class SelectiveSCAN: public NetworKit::SelectiveCommunityDetector {
public:

	SelectiveSCAN();

	virtual ~SelectiveSCAN();

	/**
	 * @param[in]	G		the graph
	 * @param[in]	seed	seed node
	 *
	 * @param[out]			the community as a set of nodes
	 */
	virtual std::unordered_set<node> run(Graph& G, node seed);
};

} /* namespace NetworKit */
#endif /* SELECTIVESCAN_H_ */
