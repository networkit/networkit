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

/**
 * Generates seed nodes.
 */
class SeedSetGenerator {

public:

	/**
	 * @param[in] G		pointer to the current graph
	 */
	SeedSetGenerator(const Graph& G);

	virtual ~SeedSetGenerator();

	/**
	 * @param[in] k 	number of required seed nodes
	 */
	virtual std::unordered_set<node> getSeeds(count k) = 0;

protected:

	const Graph& G;  	//!< pointer to current graph

};

} /* namespace NetworKit */
#endif /* SEEDSETGENERATOR_H_ */
