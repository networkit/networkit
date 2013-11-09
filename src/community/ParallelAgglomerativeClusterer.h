/*
 * ParallelAgglomerativeClusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef PARALLELAGGLOMERATIVECLUSTERER_H_
#define PARALLELAGGLOMERATIVECLUSTERER_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * A parallel agglomerative community detection algorithm, maximizing modularity.
 */
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
