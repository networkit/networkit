/*
 * ForestFireGenerator.cpp
 *
 *  Created on: 17.01.2014
 *      Author: cls
 */

#include "ForestFireGenerator.h"
#include "../auxiliary/Random.h"
#include <unordered_map>


namespace NetworKit {

ForestFireGenerator::ForestFireGenerator(double p) :p(p) {
}

std::vector<GraphEvent> ForestFireGenerator::generate(count nSteps) {

	std::unordered_map<node, bool> visited;


	auto x = [&]() {
		return Aux::Random::binomial(1 / (p * (1-p)), p);
	};

	// select k "edges" of node u
	auto selectEdges = [&](count k, node u) {
		count i = 0;
		std::set<node> edges;
		G.forNeighborsOf(u, [&](node v) {
			if (i < k) {
				edges.insert(v);
				i++;
			}
		});
		return edges;
	};

	auto selectNeighbors = [&](count k, node u) {
		std::set<node> neighbors;
		for (node c : selectEdges(k, u)) {
			if (! visited[c]) {
				neighbors.insert(c);
				visited[c] = true;
			}
		}
		return neighbors;
	};

	/** 
	 * Union of multiple sets.
	 */
	auto unite = [](std::vector<std::set<node> > sets) {
		std::set<node> all;
		for (auto set : sets) {
			all.insert(set.begin(), set.end());
		}
		return all;
	};

	std::vector<GraphEvent> stream;
	for (index step = 0; step < nSteps; ++step) {
		// select ambassador node
		node a = none;
		do {
			a = Aux::Random::integer(G.upperNodeIdBound());
		} while (! G.hasNode(a));
		assert (a != none);

		std::set<node> candidates;
		visited[a] = true;

		node s = a;
		do {
			auto neighbors = selectNeighbors(x(), a);
		} while (false); // TODO


		node v = G.addNode();
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, v));
	}
}


} /* namespace NetworKit */

