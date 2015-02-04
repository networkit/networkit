/*
 * IndependentSetFinder.cpp
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "IndependentSetFinder.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

std::string IndependentSetFinder::toString() const {
	return "TODO: implement IndependentSetFinder.toString";
}

bool IndependentSetFinder::isIndependentSet(const std::vector<bool>& set, const Graph& G) const {
	bool result = true;
	G.forEdges([&](node u, node v) {
		if (u != v) { // exclude self-loop case
			if (set[u] & set[v]) {
				DEBUG("connected nodes " , u , " and " , v , " are in the set");
				result = false;
			}
		}

	});
	return result;
}


} /* namespace NetworKit */
