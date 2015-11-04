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
#include "../matching/PathGrowingMatcher.h"
#include "../community/PLM.h"
#include "../distance/AlgebraicDistance.h"


namespace NetworKit {


MultiscaleGenerator::MultiscaleGenerator(const Graph& original) : original(original) {

}

Graph MultiscaleGenerator::generate() {



	// @param[in]	u_	coarse node
	auto replicateSubgraph = [&](node u_) {
		std::map<node, node> localNodeMap;


	};

	// V-cycle of coarsening and uncoarsening
	std::vector<Graph> coarse;
	std::vector<Graph> fine;
	std::vector<std::vector<node>> nodeMapping;
	std::vector<std::map<node, std::vector<node>>> reverseNodeMapping;

	coarse[0] = original; 	// FIXME: possibly avoid copy of the graph

	for (index level = 0; level < maxLevels; ++level) {

		// coarsen graph
		std::unique_ptr<GraphCoarsening> coarsening;

		//      select aggregation scheme
		if (aggregationScheme == "matching") {
			// TODO: select edge weighting scheme
			AlgebraicDistance ad(coarse[level]);
			ad.preprocess();
			PathGrowingMatcher matcher(coarse[level], ad.getEdgeAttribute());
			matcher.run();
			Matching matching = matcher.getMatching();
			coarsening.reset(new MatchingCoarsening(original, matching));
		} else if (aggregationScheme == "communities") {
			PLM plm(original, false, 1.0, "balanced", 32, false, false);	// recurse = false
			plm.run();
			coarsening.reset(new ParallelPartitionCoarsening(original, plm.getPartition()));
		}
		coarsening->run();


		coarse[level] = coarsening->getCoarseGraph();
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
