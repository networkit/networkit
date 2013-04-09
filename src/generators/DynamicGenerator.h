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

	GraphEventProxy* proxy;	//!< receives events produced by the generator and forwards them
	Graph* G;

public:

	DynamicGraphGenerator();

	virtual ~DynamicGraphGenerator();

	/*
	 * Set the graph proxy object.
	 */
	virtual void setProxy(GraphEventProxy& proxy);

	/*
	 * Send graph events to the proxy until termination function becomes true.
	 */
	virtual void generate(std::function<bool(void)> terminate) = 0;
};

} /* namespace NetworKit */
#endif /* DYNAMICGRAPHGENERATOR_H_ */
