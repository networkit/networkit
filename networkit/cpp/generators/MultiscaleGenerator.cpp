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


namespace NetworKit {


MultiscaleGenerator::MultiscaleGenerator(const Graph& original) : original(original) {

}

Graph MultiscaleGenerator::generate() {

	std::unique_ptr<GraphCoarsening> coarseningScheme;

	// coarsen graph
	//      aggregation scheme: return Partition
	//
	// coarse level edits
	//
	//
	// TODO:
	return Graph();
}


} /* namespace NetworKit */
