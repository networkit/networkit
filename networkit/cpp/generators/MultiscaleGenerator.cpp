/*
 * MultiscaleGenerator.h
 *
 *  Created on: Oct 27, 2015
 *      Author: Christian Staudt
 */

#include "MultiscaleGenerator.h"

#include <memory>
#include "../coarsening/GraphCoarsening.h"
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../coarsening/MatchingCoarsening.h"
#include "../matching/Matching.h"
#include "../matching/LocalMaxMatcher.h"
#include "../community/PLM.h"


namespace NetworKit {


MultiscaleGenerator::MultiscaleGenerator(const Graph& original) : original(original) {

}

Graph MultiscaleGenerator::generate() {



	// @param[in]	u_	coarse node
	auto replicateSubgraph = [&](node u_) {
		std::map<node, node> localNodeMap;


	};

	// V-cycle of coarsening and uncoarsening
	std::vector<Graph> down;
	std::vector<Graph> up;
	std::vector<std::vector<node>> nodeMapping;
	std::vector<std::map<node, std::vector<node>>> reverseNodeMapping;

	for (index level = 0; level < maxLevels; ++level) {

		// coarsen graph
		std::unique_ptr<GraphCoarsening> coarsening;

		//      select aggregation scheme
		if (aggregationScheme == "matching") {
			// TODO: select edge weighting scheme
			// PathGrowingMatcher matcher(original);
			// Matching matching = matcher.run();
			// coarsening.reset(new MatchingCoarsening(original, matching));
		} else if (aggregationScheme == "communities") {
			PLM plm(original, false, 1.0, "balanced", 32, false, false);	// recurse = false
			plm.run();
			coarsening.reset(new ParallelPartitionCoarsening(original, plm.getPartition()));
		}
		coarsening->run();


		down[level] = coarsening->getCoarseGraph();
		nodeMapping[level] = coarsening->getFineToCoarseNodeMapping();	// fine node -> coarse node
		reverseNodeMapping[level] = coarsening->getCoarseToFineNodeMapping();	//	coarse node -> collection of fine nodes

		// TODO: coarsest-level edits: delete nodes, add nodes

		// TODO: editing parameters, growth/shrink

	}




	//
	// coarse level edits
	//
	//
	// TODO: return replica
	return original;
}


} /* namespace NetworKit */
