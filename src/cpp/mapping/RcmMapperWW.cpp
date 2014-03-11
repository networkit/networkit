/*
 * RcmMapper.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#include "RcmMapperWW.h"
#include <limits>
#include <queue>

namespace NetworKit {

RcmMapperWW::RcmMapperWW() {
	// TODO Auto-generated constructor stub

}

RcmMapperWW::~RcmMapperWW() {
	// TODO Auto-generated destructor stub
}

std::map<index, index> RcmMapperWW::run(Graph& guest, Graph& host) {
	permutation pi1 = map(guest);
	permutation pi2 = map(host);

	permutation f = compose(invert(pi1), pi2);

	// create map of permutation
	std::map<index, index> mapping;
	for (index i = 0; i < f.size(); i++) {
		mapping.insert(std::pair<index, index>(i, f.at(i)));
	}

	return mapping;
}

permutation RcmMapperWW::invert(permutation pi) {
	permutation inverted_pi(pi.size());
	for (index i = 0; i < pi.size(); i++) {
		inverted_pi.at(pi.at(i)) = i;
	}

	return inverted_pi;
}

/**
 * Returns pi2 o pi1
 */
permutation RcmMapperWW::compose(permutation pi1, permutation pi2) {
	permutation composition(pi1.size());
	for (index i = 0; i < pi1.size(); i++) {
		composition.at(i) = pi2.at(pi1.at(i));
	}

	return composition;
}

// assumption: g has only one connected component
permutation RcmMapperWW::map(Graph& g) {
	// permutation
	permutation pi;
	std::queue<node> q;
	std::vector<bool> visited(g.numberOfNodes(), false);

	// find node with minimum node degree
	node s;
	edgeweight minDegree = std::numeric_limits<edgeweight>::max();
	g.forNodes([&](node v){
		if (g.weightedDegree(v) < minDegree) {
			s = v;
			minDegree = g.weightedDegree(v);
		}
	});


	// insert node with minimum degree
	q.push(s);
	visited.at(s) = true;

	while (!q.empty()) {
		node v = q.front();
		q.pop();
		pi.push_back(v);

		// visit neighbours of v in degree increasing order
		Comparator c(g);
		std::priority_queue<node, std::vector<node>, Comparator> neighbours(c);
		g.forNeighborsOf(v, [&](node u) {
			if (visited.at(u)) return;
			neighbours.push(u);
		});

		while (!neighbours.empty()) {
			q.push(neighbours.top());
			visited.at(neighbours.top()) = true;
			neighbours.pop();
		}
	}

//	std::reverse(pi.begin(), pi.end());

	return pi;
}


} /* namespace NetworKit */
