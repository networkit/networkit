/*
 * RandomWalkSeedSet.cpp
 *
 *  Created on: 20.06.2013
 *      Author: cls
 */

#include "RandomWalkSeedSet.h"

namespace NetworKit {

RandomWalkSeedSet::RandomWalkSeedSet(const Graph& G, count nSteps) : SeedSetGenerator(G), nSteps(nSteps) {
	// TODO Auto-generated constructor stub

}

RandomWalkSeedSet::~RandomWalkSeedSet() {
	// TODO Auto-generated destructor stub
}

std::unordered_set<node> RandomWalkSeedSet::getSeeds(count k) {
	std::unordered_set<node> S;


	auto randomWalkStep = [&](node previous) {
		std::vector<node> neighbors;
		G.forNeighborsOf(previous, [&](node v){
			neighbors.push_back(v);
		});
		if (neighbors.size() > 0) {
			index ri = randInt.generate(0, neighbors.size() - 1);
			node next = neighbors.at(ri);
			return next;
		} else {
			throw std::runtime_error("cannot perform random walk on isolated node");
		}
	};

	// find random initial node
	index z = G.upperNodeIdBound();

	node r;
	do {
		r = randInt.generate(0, z);
	} while (! G.hasNode(r) && (G.degree(r) > 0));

	S.insert(r);

	node current = r;
	while (S.size() < k) {
		// perform random walk step
		for (count i = 0; i < nSteps; ++i){
			current = randomWalkStep(current);
		}
		S.insert(current);
	}

	return S;
}

} /* namespace NetworKit */
