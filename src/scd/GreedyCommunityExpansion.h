/*
 * GreedyCommunityExpansion.h

 *
 *  Created on: 03.06.2013
 *      Author: cls , Yassine Marrakchi
 */

#ifndef GREEDYCOMMUNITYEXPANSION_H_
#define GREEDYCOMMUNITYEXPANSION_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"
#include "CommunityTrimming.h"
#include "QualityObjective.h"
#include "Acceptability.h"


namespace NetworKit {


/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
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

	virtual std::map<node, std::unordered_set<node>> seedSetExpansion(Graph& G, std::vector<node> set);

	virtual double clusterClusterSimilarity (std::unordered_set<node>& community1, std::unordered_set<node>& community2);
};

} /* namespace NetworKit */
#endif /* GREEDYCOMMUNITYEXPANSION_H_ */



