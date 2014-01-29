/*
 * ForestFireGenerator.cpp
 *
 *  Created on: 17.01.2014
 *      Author: cls
 */

#include "ForestFireGenerator.h"
#include "../auxiliary/Random.h"


namespace NetworKit {

ForestFireGenerator::ForestFireGenerator(double p) :p(p) {
}

std::vector<GraphEvent> ForestFireGenerator::generate(count nSteps) {

	auto x = [&]() {
		return Aux::Random::binomial(1 / (p * (1-p)), p);
	};

	// select k neighbors of node u
	auto selectNeighbors = [&](count k, node u) {
		count i = 0;
		std::set<node> neighbors;
		G.forNeighborsOf(u, [&](node v) {
			if (i < k) {
				neighbors.insert(v);
				i++;
			}
		});
		return neighbors;
	};

	std::vector<GraphEvent> stream;
	for (index step = 0; step < nSteps; ++step) {
		// select ambassador node
		node w = Aux::Random::integer(G.upperNodeIdBound());

		std::set<node> candidates;

		for (node u : selectNeighbors(x(), w)) {
			// TODO:
		}


		node v = G.addNode();
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, v));
	}
}


} /* namespace NetworKit */

