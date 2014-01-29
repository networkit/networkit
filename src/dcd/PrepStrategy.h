/*
 * PrepStrategy.h
 *
 *  Created on: 18.04.2013
 *      Author: cls
 */

#ifndef PREPSTRATEGY_H_
#define PREPSTRATEGY_H_

#include "../dynamics/GraphEventHandler.h"

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

	virtual void onTimeStep();

	/**
	 * String representation - PrepStrategy subclasses must implement this.
	 */
	virtual std::string toString() = 0;

};

} /* namespace NetworKit */
#endif /* PREPSTRATEGY_H_ */
