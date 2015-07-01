/*
 * Assortativity.cpp
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#include "Assortativity.h"

namespace NetworKit {

Assortativity::Assortativity(const Graph& G, const std::vector<double>& attribute) : G(G), attribute(attribute) {
	if (attribute.size() < G.upperNodeIdBound()) {
		throw std::runtime_error("attribute list has incorrect length: there must be an entry for each node");
	}
}


void Assortativity::run() {
	throw std::runtime_error("TODO");
}


double Assortativity::getCoefficient() const {
	return 0.0;
}


}
