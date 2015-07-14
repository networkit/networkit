/*
 * GraphEventHandler.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTHANDLER_H_
#define GRAPHEVENTHANDLER_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class GraphEventHandler {

public:
	virtual void onNodeAddition(node u) = 0;

	virtual void onNodeRemoval(node u) = 0;

	virtual void onNodeRestoration(node u) = 0;

	virtual void onEdgeAddition(node u, node v, edgeweight w = 1.0) = 0;

	virtual void onEdgeRemoval(node u, node v, edgeweight w = 1.0) = 0;

	virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) = 0;
	
	virtual void onWeightIncrement(node u, node v, edgeweight wOld, edgeweight delta) = 0;

	virtual void onTimeStep() = 0;
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTHANDLER_H_ */
