/*
 * PrepStrategy.h
 *
 *  Created on: 18.04.2013
 *      Author: cls
 */

#ifndef PREPSTRATEGY_H_
#define PREPSTRATEGY_H_

#include "../dynamics/GraphEventHandler.h"
#include "DynamicCommunityDetector.h"

namespace NetworKit {

/**
 * Abstract base class for all dynamic community detection prep strategies.
 * A prep strategy handles incoming graph events, and e.g. producing a preclustering
 * for the community detection algorithm.
 *
 */
class PrepStrategy: public NetworKit::GraphEventHandler {


public:

	PrepStrategy();

	virtual ~PrepStrategy();

	virtual void onNodeAddition(node u) = 0;

	virtual void onNodeRemoval(node u) = 0;

	virtual void onEdgeAddition(node u, node v) = 0;

	virtual void onEdgeRemoval(node u, node v) = 0;

	virtual void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) = 0;

};

} /* namespace NetworKit */
#endif /* PREPSTRATEGY_H_ */
