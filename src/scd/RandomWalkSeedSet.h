/*
 * RandomWalkSeedSet.h
 *
 *  Created on: 20.06.2013
 *      Author: cls
 */

#ifndef RANDOMWALKSEEDSET_H_
#define RANDOMWALKSEEDSET_H_

#include "SeedSetGenerator.h"

#include "../auxiliary/RandomInteger.h"
#include "../auxiliary/RandomProbability.h"


namespace NetworKit {

class RandomWalkSeedSet: public NetworKit::SeedSetGenerator {
public:

	RandomWalkSeedSet(const Graph& G, count nSteps = 2);

	virtual ~RandomWalkSeedSet();

	virtual std::unordered_set<node> getSeeds(count k);

protected:

	count nSteps;
	Aux::RandomInteger randInt;
	Aux::RandomProbability randProb;

};

} /* namespace NetworKit */
#endif /* RANDOMWALKSEEDSET_H_ */
