/*
 * ParallelAgglomerativeClusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef PARALLELAGGLOMERATIVECLUSTERER_H_
#define PARALLELAGGLOMERATIVECLUSTERER_H_

#include "Clusterer.h"
#include "../scoring/ModularityScoring.h"
#include "../matching/ParallelMatcher.h"
#include "../matching/PathGrowingMatcher.h"
#include "../coarsening/MatchingContracter.h"
#include "../coarsening/ClusteringProjector.h"
#include "../Globals.h"


namespace NetworKit {

class ParallelAgglomerativeClusterer: public NetworKit::Clusterer {

public:

	ParallelAgglomerativeClusterer();

	virtual ~ParallelAgglomerativeClusterer();

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;

	virtual Clustering run(Graph& G);
};

} /* namespace NetworKit */
#endif /* SCOREMATCHCONTRACT_H_ */
