/*
 * RandomClusterer.h
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef RANDOMCLUSTERER_H_
#define RANDOMCLUSTERER_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * Dummy community detection algorithm: Creates a random partition of a network.
 */
class RandomClusterer: public NetworKit::Clusterer {
public:
	RandomClusterer();
	virtual ~RandomClusterer();

	virtual Partition run(Graph& G);
};

} /* namespace NetworKit */
#endif /* RANDOMCLUSTERER_H_ */
