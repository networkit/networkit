/*
 * DynamicGenerator.h
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#ifndef DYNAMICGRAPHGENERATOR_H_
#define DYNAMICGRAPHGENERATOR_H_

#include <functional>

#include "../dynamics/GraphEventProxy.h"

namespace NetworKit {

class DynamicGraphGenerator {

protected:

	GraphEventProxy* Gproxy;	//!< receives events produced by the generator and forwards them
	Graph* G;

public:

	DynamicGraphGenerator(GraphEventProxy& proxy);

	virtual ~DynamicGraphGenerator();

	/**
	 * The generator may expect the graph to be in a certain initial state. Call this method first.
	 */
	virtual void initializeGraph() = 0;

	/*
	 * Send graph events to the proxy while function does not return false.
	 */
	virtual void generateWhile(std::function<bool(void)> cont) = 0;
};

} /* namespace NetworKit */
#endif /* DYNAMICGRAPHGENERATOR_H_ */
