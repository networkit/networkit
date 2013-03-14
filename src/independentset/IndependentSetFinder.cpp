/*
 * IndependentSetFinder.cpp
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "IndependentSetFinder.h"

namespace EnsembleClustering {

IndependentSetFinder::IndependentSetFinder() {
	// TODO Auto-generated constructor stub

}

IndependentSetFinder::~IndependentSetFinder() {
	// TODO Auto-generated destructor stub
}

std::string IndependentSetFinder::toString() const {
	return "TODO: implement IndependentSetFinder.toString";
}

bool IndependentSetFinder::isIndependentSet(const std::vector<bool>& set, const Graph& G) const {
	G.parallelForEdges([&](node u, node v) {
		if (set[u] & set[v]) {
			return false;
		}
	});
	return true;
}


} /* namespace EnsembleClustering */
