/*
 * SeedSetGenerator.h
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef SEEDSETGENERATOR_H_
#define SEEDSETGENERATOR_H_

#include "../graph/Graph.h"

namespace NetworKit {

class SeedSetGenerator {

public:

	SeedSetGenerator(Graph& G);

	virtual ~SeedSetGenerator();

	virtual std::vector<node> getSeeds() = 0;

protected:

	Graph* G;

};

} /* namespace NetworKit */
#endif /* SEEDSETGENERATOR_H_ */
