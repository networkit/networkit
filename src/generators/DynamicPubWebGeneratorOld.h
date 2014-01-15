/*
 * DynamicPubWebGenerator.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef DYNAMICPUBWEBGENERATOR_OLD_H_
#define DYNAMICPUBWEBGENERATOR_OLD_H_

#include "DynamicGraphSource.h"
#include "PubWebGenerator.h"

namespace NetworKit {

class DynamicPubWebGeneratorOld: public DynamicGraphSource {
protected:
//	PubWebGenerator staticGen;

	void determineNeighbors(Graph& g);
	void determineNeighborsOf(Graph& g, node u);
	void moveNodeIntoUnitSquare(float& x, float& y);
	float squaredDistanceInUnitTorus(float x1, float y1, float x2, float y2);
	void chooseDenseAreaSizes();
	void fillDenseAreas(Graph& g);
	void spreadRemainingNodes(Graph& g);
	void chooseClusterSizes();
	void addNodesToArea(index area, count num, Graph& g);
	bool isValidEdge(Graph& g, node u, node v);

	void moveNodesRandomly();

public:
	DynamicPubWebGeneratorOld(count numInitialNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors);
	virtual ~DynamicPubWebGeneratorOld();

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
