/*
 * RandomSeedSet.cpp
 *
 *  Created on: 20.06.2013
 *      Author: cls
 */

#include "RandomSeedSet.h"

namespace NetworKit {

RandomSeedSet::RandomSeedSet(const Graph& G) : SeedSetGenerator(G) {
	// TODO Auto-generated constructor stub

}

RandomSeedSet::~RandomSeedSet() {
	// TODO Auto-generated destructor stub
}

std::unordered_set<node> RandomSeedSet::getSeeds(count k) {
	if (G.numberOfNodes() < k) {
		throw std::runtime_error("Graph has fewer nodes than requested seed set");
	}

	std::unordered_set<node> S;
	while (S.size() < k) {
		node r = Aux::RandomInteger::generate(0, G.upperNodeIdBound() - 1);
		if (G.hasNode(r)) {
			S.insert(r);
		}
	}

	return S;

}

} /* namespace NetworKit */
