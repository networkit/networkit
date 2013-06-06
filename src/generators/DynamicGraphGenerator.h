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


	/**
	 * Perform one generative step - as defined by the implementation.
	 */
	virtual void generate() = 0;

	/*
	 * Continue generating while function does not return false.
	 * @param[in]	cont	generator continues when this function returns true
	 */
	virtual void generateWhile(std::function<bool(void)> cont);

	/**
	 * Continue generating until the number of nodes reaches this upper limit.
	 * @param[in]	n	number of nodes
	 */
	virtual void generateNodes(count n);


	/**
	 * Continue generating until the number of edges reaches this upper limit.
	 * @param[in]	m	number of edges
	 */
	virtual void generateEdges(count m);
};

} /* namespace NetworKit */
#endif /* DYNAMICGRAPHGENERATOR_H_ */
