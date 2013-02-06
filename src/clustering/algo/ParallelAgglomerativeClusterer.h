/*
 * ParallelAgglomerativeClusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef PARALLELAGGLOMERATIVECLUSTERER_H_
#define PARALLELAGGLOMERATIVECLUSTERER_H_

#include "Clusterer.h"
#include "../../scoring/ModularityScoring.h"
#include "../../matching/ParallelMatcher.h"
#include "../../coarsening/MatchingContracter.h"


namespace EnsembleClustering {

class ParallelAgglomerativeClusterer: public EnsembleClustering::Clusterer {

public:

	ParallelAgglomerativeClusterer();

	virtual ~ParallelAgglomerativeClusterer();

	virtual Clustering run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* SCOREMATCHCONTRACT_H_ */
