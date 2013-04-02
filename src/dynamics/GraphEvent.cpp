/*
 * GraphEvent.cpp
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#include "GraphEvent.h"

namespace NetworKit {

GraphEvent::GraphEvent(GraphEventType type, node u, node v, edgeweight w) : type(type), u(u), v(v), w(w) {
}

GraphEvent::~GraphEvent() {
	// TODO Auto-generated destructor stub
}

} /* namespace NetworKit */
