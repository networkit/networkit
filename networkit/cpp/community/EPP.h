/*
 * EnsemblePreprocessing.h
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ENSEMBLEPREPROCESSING_H_
#define ENSEMBLEPREPROCESSING_H_

#include <vector>
#include <memory>
#include "../community/CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"
#include "../overlap/Overlapper.h"


namespace NetworKit {

/**
 * @ingroup community
 * EPP - Ensemble Preprocessing community detection algorithm.
 * Combines multiple base algorithms and a final algorithm. A consensus of the
 * solutions of the base algorithms is formed and the graph is coarsened accordingly.
 * Then the final algorithm operates on the coarse graph and determines a solution
 * for the input graph.
 */
class EPP: public NetworKit::CommunityDetectionAlgorithm {

protected:

	std::unique_ptr<CommunityDetectionAlgorithm> finalClusterer;	//!< final clustering algorithm
	std::vector<std::unique_ptr<CommunityDetectionAlgorithm>> baseClusterers; //!< ensemble of base clusterers

	std::unique_ptr<Overlapper> overlap; //!< clustering overlap algorithm

	Partition core;
	std::vector<Partition> baseClusterings;

public:
	/**
	 * Constructor to the EPP community detection algorithm.
	 *
	 * @param[in]	G	input graph
	 */
	EPP(const Graph& G);

	/**
	 * Add a base clusterer to the ensemble.
	 *
	 * @param base A base clusterer.
	 */
	virtual void addBaseClusterer(std::unique_ptr<CommunityDetectionAlgorithm>& base);


	/**
	 * Set final clusterer.
	 *
	 * @param final The final clusterer.
	 */
	virtual void setFinalClusterer(std::unique_ptr<CommunityDetectionAlgorithm>& final);

	/**
	 * Set overlap algorithm which combines the results of the base clusterers.
	 *
	 * @param overlap The overlap algorithm.
	 */
	virtual void setOverlapper(std::unique_ptr<Overlapper>& overlap);

	/**
	 * Run the ensemble clusterer.
	 */
	virtual void run();

	/**
	 * String representation of EPP class.
	 * @return string representation.
	 */
	virtual std::string toString() const;


	std::vector<Partition> getBasePartitions() const;

	Partition getCorePartition() const;

};

} /* namespace NetworKit */
#endif /* ENSEMBLEPREPROCESSING_H_ */
