/*
 * RandomSeedSet.h
 *
 *  Created on: 20.06.2013
 *      Author: cls
 */

#ifndef RANDOMSEEDSET_H_
#define RANDOMSEEDSET_H_

#include "SeedSetGenerator.h"

#include "../auxiliary/RandomInteger.h"

namespace NetworKit {

class RandomSeedSet: public NetworKit::SeedSetGenerator {

public:

	RandomSeedSet(const Graph& G);

	virtual ~RandomSeedSet();

	/**
	 * Get k random nodes from the graph.
	 */
	virtual std::unordered_set<node> getSeeds(count k);

};

} /* namespace NetworKit */
#endif /* RANDOMSEEDSET_H_ */
