/*
 * GraphEventSource.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTSOURCE_H_
#define GRAPHEVENTSOURCE_H_

namespace NetworKit {

#include "GraphEvent.h"

class GraphEventSource {

public:

	GraphEventSource();

	virtual ~GraphEventSource();

	virtual GraphEvent emit() = 0;
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTSOURCE_H_ */
