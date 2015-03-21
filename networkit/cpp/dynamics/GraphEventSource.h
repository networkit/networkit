/*
 * GraphEventSource.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTSOURCE_H_
#define GRAPHEVENTSOURCE_H_

#include "GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class GraphEventSource {

public:
	virtual GraphEvent emit() = 0;
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTSOURCE_H_ */
