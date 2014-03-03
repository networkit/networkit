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

ForestFireGenerator::ForestFireGenerator(double p) :p(p), firstCall(true) {
}

std::vector<GraphEvent> ForestFireGenerator::generate(count nSteps) {

	std::vector<GraphEvent> stream;
	std::set<node> empty;


	/** 
	 * Random binomially distributed number depending on "burning probability".
	 */
	auto x = [&]() {
		count k = Aux::Random::binomial(1 / (p * (1-p)), p);
		TRACE("k = ", k);
		return k;
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





	auto connectNewNode = [&]() {

		std::unordered_map<node, bool> visited;

		// select k "edges" of node u
		auto selectEdges = [&](node u, count k) {
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

		/** 
		 * Given a node, select a set of neighbors according to 
		 */
		auto select = [&](node u, count k) {
			std::set<node> neighbors;
			for (node c : selectEdges(k, u)) {
				if (! visited[c]) {
					neighbors.insert(c);
					visited[c] = true;
				}
			}
			return neighbors;
		};


		std::function<std::set<node>(std::set<node>)> foo = [&](std::set<node> V) {
			std::vector<std::set<node> > S;
			for (node v : V) {
				S.push_back(select(v, x()));
			}
			std::set<node> U = unite(S);
			if (U.empty()) {
				return empty;
			} else {
				std::set<node> Z = foo(U);
				U.insert(Z.begin(), Z.end());
				return U;
			}
		};


		// select ambassador node
		node a = none;
		do {
			a = Aux::Random::integer(G.upperNodeIdBound());
		} while (! G.hasNode(a));
		assert (a != none);
		DEBUG("selected ambassador: ", a);

		node v = G.addNode();
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, v));
		DEBUG("created node ", v);

		visited[a] = true;
		std::set<node> candidates = {a};
		std::set<node> additional = foo(candidates);
		candidates.insert(additional.begin(), additional.end());
		DEBUG("candidates: ", candidates);

		for (node c : candidates) {
			G.addEdge(v, c);
			stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v, c));
		}
	};

	// initial graph
	if (firstCall) {
		node s = G.addNode();
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, s));
		node t = G.addNode();
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, t));
		G.addEdge(s, t);
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, s, t));
		firstCall = false;
	}


	for (index step = 0; step < nSteps; ++step) {
		connectNewNode();
	}
	return stream;
}


} /* namespace NetworKit */

