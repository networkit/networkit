/*
 * GraphEventHandler.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTHANDLER_H_
#define GRAPHEVENTHANDLER_H_

namespace NetworKit {

#include "Graph.h"

class GraphEventHandler {

	void onNodeAddition(node u) = 0;

	void onNodeRemoval(node u) = 0;

	void onEdgeAddition(node u, node v) = 0;

	void onEdgeRemoval(node u, node v) = 0;

	void onWeightUpdate(node u, node v, edgeweight w);
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTHANDLER_H_ */
