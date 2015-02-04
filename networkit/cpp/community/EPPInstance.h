/*
 * EPPInstance.h
 *
 *  Created on: 9.11.2014
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef EPPINSTANCE_H_
#define EPPINSTANCE_H_

#include <vector>
#include <memory>
#include "../community/CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"
#include "../overlap/Overlapper.h"
#include "../overlap/HashingOverlapper.h"
#include "../community/PLP.h"
#include "../community/PLM.h"


namespace NetworKit {

/**
 * @ingroup community
 * EPPInstance - Ensemble Preprocessing community detection algorithm.
 * Combines multiple base algorithms and a final algorithm. A consensus of the
 * solutions of the base algorithms is formed and the graph is coarsened accordingly.
 * Then the final algorithm operates on the coarse graph and determines a solution
 * for the input graph.
 */
class EPPInstance: public NetworKit::CommunityDetectionAlgorithm {

protected:


	Partition core;
	std::vector<Partition> baseClusterings;

	count ensembleSize;

public:
	/**
	 * Constructor to the EPPInstance community detection algorithm.
	 *
	 * @param[in]	G	input graph
	 */
	EPPInstance(const Graph& G, count ensembleSize = 4);

	/**
	 * Run the ensemble clusterer.
	 */
	virtual void run();

	/**
	 * String representation of EPPInstance class.
	 * @return string representation.
	 */
	virtual std::string toString() const;


	std::vector<Partition> getBasePartitions() const;

	Partition getCorePartition() const;

};

} /* namespace NetworKit */
#endif /* ENSEMBLEPREPROCESSING_H_ */
