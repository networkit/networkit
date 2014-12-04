/*
 * GraphEvent.cpp
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#include "GraphEvent.h"

#include <sstream>

namespace NetworKit {

GraphEvent::GraphEvent(GraphEvent::Type type, node u, node v, edgeweight w) : type(type), u(u), v(v), w(w) {
}

std::string GraphEvent::toString() {
	std::stringstream ss;
	if (this->type == GraphEvent::NODE_ADDITION) {
		ss << "an(" << u << ")";
	} else if (this->type == GraphEvent::NODE_REMOVAL) {
		ss << "dn(" << u << ")";
	} else if (this->type == GraphEvent::EDGE_ADDITION) {
		ss << "ae(" << u << "," << v << "," << w << ")";
	} else if (this->type == GraphEvent::EDGE_REMOVAL) {
		ss << "de(" << u << "," << v << ")";
	} else if (this->type == GraphEvent::EDGE_WEIGHT_UPDATE) {
		ss << "ce(" << u << "," << v << ")";
	} else if (this->type == GraphEvent::TIME_STEP) {
		ss << "st";
	}
	return ss.str();
}


} /* namespace NetworKit */

