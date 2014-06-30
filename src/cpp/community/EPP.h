/*
 * EnsemblePreprocessing.h
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ENSEMBLEPREPROCESSING_H_
#define ENSEMBLEPREPROCESSING_H_

#include <vector>
#include "../community/CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"
#include "../overlap/Overlapper.h"


namespace NetworKit {

/**
 * EPP - Ensemble Preprocessing community detection algorithm.
 * Combines multiple base algorithms and a final algorithm. A consensus of the
 * solutions of the base algorithms is formed and the graph is coarsened accordingly.
 * Then the final algorithm operates on the coarse graph and determines a solution
 * for the input graph.
 */
class EPP: public NetworKit::CommunityDetectionAlgorithm {

protected:

	CommunityDetectionAlgorithm* finalClusterer;	//!< final clustering algorithm
	std::vector<CommunityDetectionAlgorithm*> baseClusterers; //!< ensemble of base clusterers

	Overlapper* overlap; //!< clustering overlap algorithm

public:

	/** Default constructor */
	EPP();

	/** Default destructor */
	virtual ~EPP() = default;

	/**
	 * Add a base clusterer to the ensemble.
	 *
	 * @param base A base clusterer.
	 */
	virtual void addBaseClusterer(CommunityDetectionAlgorithm&  base);


	/**
	 * Set final clusterer.
	 *
	 * @param final The final clusterer.
	 */
	virtual void setFinalClusterer(CommunityDetectionAlgorithm& final);

	/**
	 * Set overlap algorithm which combines the results of the base clusterers.
	 *
	 * @param overlap The overlap algorithm.
	 */
	virtual void setOverlapper(Overlapper& overlap);

	/**
	 * Run the ensemble clusterer on @a G and return the result in a Partition.
	 *
	 * @param G The graph.
	 * @return A Partition of the clustering.
	 */
	virtual Partition run(Graph& G);

	/**
	 * String representation of EPP class.
	 * @return string representation.
	 */
	virtual std::string toString() const;

};

} /* namespace NetworKit */
#endif /* ENSEMBLEPREPROCESSING_H_ */
