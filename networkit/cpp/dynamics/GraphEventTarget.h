/*
 * GraphEventTARGET.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTTARGET_H_
#define GRAPHEVENTTARGET_H_

#include "GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class GraphEventTarget {

public:

	/** Default destructor */
	virtual ~GraphEventTarget();

	virtual void receive(GraphEvent event) = 0;
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTTARGET_H_ */
