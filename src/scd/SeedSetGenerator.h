/*
 * SeedSetGenerator.h
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef SEEDSETGENERATOR_H_
#define SEEDSETGENERATOR_H_

#include <unordered_set>

#include "../graph/Graph.h"

namespace NetworKit {

class SeedSetGenerator {

public:

	SeedSetGenerator(const Graph& G);

	virtual ~SeedSetGenerator();

	virtual std::unordered_set<node> getSeeds(count k) = 0;

protected:

	const Graph& G;

};

} /* namespace NetworKit */
#endif /* SEEDSETGENERATOR_H_ */
