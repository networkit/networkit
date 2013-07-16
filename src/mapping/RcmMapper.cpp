/*
 * RcmMapper.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#include "RcmMapper.h"

namespace NetworKit {

RcmMapper::RcmMapper() {
	// TODO Auto-generated constructor stub

}

RcmMapper::~RcmMapper() {
	// TODO Auto-generated destructor stub
}

Permutation RcmMapper::invert(const Permutation& piIn) const {
	count n = piIn.size();
	Permutation piOut(n);

	for (index i = 0; i < n; ++i) {
		piOut[piIn[i]] = i;
	}

	return piOut;
}

Mapping RcmMapper::run(Graph& guest, Graph& host) {
	std::map<index, index> mapping;
	count n = guest.numberOfNodes();
	assert(host.numberOfNodes() == n); // assumption: guest and host are of equal size

	// permute guest graph according to RCM
	Permutation piGuest = permute(guest);
	Permutation piGuestInverted = invert(piGuest);

	// permute host graph according to RCM
	Permutation piHost = permute(host);

	// TODO: fill in mapping according to permutations
	for (index i = 0; i < n; ++i) {
		mapping.insert(std::make_pair(i, piHost[piGuestInverted[i]]));
	}

	return mapping;
}

Permutation RcmMapper::permute(const Graph& graph) const {
	Permutation permutation;
	std::set<index> unvisited;

	// init unvisited => insert all nodes
	graph.forNodes([&](node v) {
		unvisited.insert(v);
	});

	while (! unvisited.empty()) {
		// for each connected component
		index argmin = graph.argminDegree();
		std::queue<index> q;
		q.push(argmin);

		while (! q.empty()) {
			// process current front node
			node v = q.front();
			q.pop();
			unvisited.erase(v);
			permutation.push_back(v);

			// enqueue unvisited neighbors in degree-increasing order
			graph.forEdgesOfInDegreeIncreasingOrder(v, [&](node v, node u) {
				if (unvisited.count(u) > 0) {
					q.push(u);
				}
			});
		}
	}

	// reverse since it is RCM
	std::reverse(permutation.begin(), permutation.end());

	return permutation;
}

} /* namespace NetworKit */
