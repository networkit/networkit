/*
 * TNodeDistance.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "TNodeDistance.h"

namespace NetworKit {

TNodeDistance::TNodeDistance(const Graph& G) : G(G) {
	// TODO Auto-generated constructor stub

}

TNodeDistance::~TNodeDistance() {
	// TODO Auto-generated destructor stub
}

void TNodeDistance::initialize(const Parameters& param) {
	// do nothing
}

double TNodeDistance::distance(node u, node v) {
	throw std::runtime_error("abstract method called - this needs to be overwritten in subclass");
}

} /* namespace NetworKit */
