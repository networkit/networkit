/*
 * StaticGraphGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef STATICGRAPHGENERATOR_H_
#define STATICGRAPHGENERATOR_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Abstract base class for static graph generators.
 */
class StaticGraphGenerator {

public:

	StaticGraphGenerator();

	virtual ~StaticGraphGenerator();

	virtual Graph generate() = 0;

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _generate() {
		return new Graph{std::move(generate())};
	};
};

} /* namespace NetworKit */
#endif /* STATICGRAPHGENERATOR_H_ */
