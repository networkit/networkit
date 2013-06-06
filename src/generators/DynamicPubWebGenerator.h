/*
 * DynamicPubWebGenerator.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef DYNAMICPUBWEBGENERATOR_H_
#define DYNAMICPUBWEBGENERATOR_H_

#include "DynamicGraphGenerator.h"
#include "PubWebGenerator.h"

namespace NetworKit {

class DynamicPubWebGenerator: public DynamicGraphGenerator {
protected:
	PubWebGenerator staticGen;

public:
	DynamicPubWebGenerator(GraphEventProxy& proxy, count numInitialNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);
	virtual ~DynamicPubWebGenerator();

	/**
	 * The generator may expect the graph to be in a certain initial state. Call this method first.
	 */
	virtual void initializeGraph();

	/*
	 * Send graph events to the proxy until termination function becomes true.
	 */
	virtual void generateWhile(std::function<bool(void)> cont);

	node addNode();
};

} /* namespace NetworKit */
#endif /* DYNAMICPUBWEBGENERATOR_H_ */
