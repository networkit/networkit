/*
 * EnsemblePreprocessing.h
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef ENSEMBLEPREPROCESSING_H_
#define ENSEMBLEPREPROCESSING_H_

#include <vector>
#include "../community/Clusterer.h"
#include "../clustering/Clustering.h"
#include "../overlap/Overlapper.h"


namespace NetworKit {

/**
 * EPP - Ensemble Preprocessing community detection algorithm.
 * Combines multiple base algorithms and a final algorithm. A consensus of the
 * solutions of the base algorithms is formed and the graph is coarsened accordingly.
 * Then the final algorithm operates on the coarse graph and determines a solution
 * for the input graph.
 */
class EPP: public NetworKit::Clusterer {

protected:

	Clusterer* finalClusterer;	//!< final clustering algorithm
	std::vector<Clusterer*> baseClusterers; //!< ensemble of base clusterers

	Overlapper* overlap; //!< clustering overlap algorithm

public:

	EPP();

	virtual ~EPP() = default;

	/**
	 * Add a base clusterer to the ensemble.
	 */
	virtual void addBaseClusterer(Clusterer&  base);


	/**
	 * Set final clusterer.
	 */
	virtual void setFinalClusterer(Clusterer& final);

	/**
	 * Set overlap algorithm which combines the results of the base clusterers.
	 */
	virtual void setOverlapper(Overlapper& overlap);

	/**
	 * Run the ensemble clusterer.
	 */
	virtual Clustering run(Graph& G);

	/**
	 * @return string representation.
	 */
	virtual std::string toString() const;

};

} /* namespace NetworKit */
#endif /* ENSEMBLEPREPROCESSING_H_ */
