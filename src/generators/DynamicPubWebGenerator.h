/*
 * DynamicPubWebGenerator.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef DYNAMICPUBWEBGENERATOR_H_
#define DYNAMICPUBWEBGENERATOR_H_

#include "DynamicGraphSource.h"
#include "PubWebGenerator.h"

namespace NetworKit {

class DynamicPubWebGenerator: public DynamicGraphSource {
protected:
	PubWebGenerator staticGen;

	void moveNodesRandomly();

public:
	DynamicPubWebGenerator(count numInitialNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);
	virtual ~DynamicPubWebGenerator();

	/**
	 * The generator may expect the graph to be in a certain initial state. Call this method first.
	 */
	virtual void initializeGraph();

	/*
	 * Send graph events to the proxy until termination function becomes true.
	 */
	virtual void generate();

	node addNode();
};

} /* namespace NetworKit */
#endif /* DYNAMICPUBWEBGENERATOR_H_ */
