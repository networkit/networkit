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
 * @ingroup generators
 * Abstract base class for static graph generators.
 */
class StaticGraphGenerator {

public:

	/** Default destructor */
	virtual ~StaticGraphGenerator() = default;

	virtual Graph generate() = 0;
};

} /* namespace NetworKit */
#endif /* STATICGRAPHGENERATOR_H_ */
