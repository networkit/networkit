/*
 * GreedyCommunityExpansion.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GREEDYCOMMUNITYEXPANSION_H_
#define GREEDYCOMMUNITYEXPANSION_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"


namespace NetworKit {

class GreedyCommunityExpansion: public NetworKit::SelectiveCommunityDetector {

public:

	GreedyCommunityExpansion();

	virtual ~GreedyCommunityExpansion();

	/**
	 * @param[in]	G	input graph
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	virtual std::unordered_set<node> run(Graph& G, node s);
};

} /* namespace NetworKit */
#endif /* GREEDYCOMMUNITYEXPANSION_H_ */
