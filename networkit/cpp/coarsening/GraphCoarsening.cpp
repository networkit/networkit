/*
 * GraphCoarsening.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "GraphCoarsening.h"

namespace NetworKit {

GraphCoarsening::GraphCoarsening(const Graph& G) : Algorithm(), G(G) {

}

Graph GraphCoarsening::getCoarseGraph() const {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return Gcoarsed;
}

std::vector<node> GraphCoarsening::getNodeMapping() const {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return nodeMapping;
}

std::string GraphCoarsening::toString() const {
	return "GraphCoarsening base class";
}

}
