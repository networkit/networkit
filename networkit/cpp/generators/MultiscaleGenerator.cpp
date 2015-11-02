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
#include "../community/PLP.h"


namespace NetworKit {


MultiscaleGenerator::MultiscaleGenerator(const Graph& original) : original(original) {

}

Graph MultiscaleGenerator::generate() {

	std::vector<Graph> coarseGraphs;
	std::vector<Graph> fineGraphs;

	// coarsen graph
	//      aggregation scheme
	std::unique_ptr<GraphCoarsening> coarsening;

	if (aggregationScheme == "matching") {
		// LocalMaxMatcher matcher(original);
		// Matching matching = matcher.run();
		// coarsening.reset(new MatchingCoarsening(original, matching));
	} else if (aggregationScheme == "communities") {
		PLP plp(original);
		plp.run();
		coarsening.reset(new ParallelPartitionCoarsening(original, plp.getPartition()));
	}
	coarsening->run();


	Graph C = coarsening->getCoarseGraph();


	//
	// coarse level edits
	//
	//
	// TODO:
	return C;
}


} /* namespace NetworKit */
