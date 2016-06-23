/*
 * DynamicGraphGenerator.h
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef DYNAMICGRAPHGENERATOR_H_
#define DYNAMICGRAPHGENERATOR_H_

#include "../graph/Graph.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup generators
 * Abstract base class for a dynamic graph generator (in the new dynamic architecture).
 * The generator produces a stream of events.
 */
class DynamicGraphGenerator {

public:
	/** Default destructor */
	virtual ~DynamicGraphGenerator() = default;

	/**
	 * Generate event stream.
	 *
	 * @param[in]	nSteps	number of time steps in the event stream
	 */
	virtual std::vector<GraphEvent> generate(count nSteps) = 0;

protected:

	Graph G; // the graph instance
};

} /* namespace NetworKit */

#endif /* DYNAMICGRAPHGENERATOR_H_ */
