/*
 * RcmMapper.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#include "RcmMapper.h"

namespace NetworKit {

RcmMapper::RcmMapper() {

}

RcmMapper::~RcmMapper() {

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
	Mapping mapping;
	count n = guest.numberOfNodes();
	assert(host.numberOfNodes() == n); // assumption: guest and host are of equal size

	// permute guest graph according to RCM
	Permutation piGuest = permute(guest);
	Permutation piGuestInverted = invert(piGuest);

	// permute host graph according to RCM
	Permutation piHost = permute(host);

	// fill in mapping according to permutations
	for (index i = 0; i < n; ++i) {
		mapping.insert(std::make_pair(i, piHost[piGuestInverted[i]]));
	}

	return mapping;
}

Permutation RcmMapper::permute(const Graph& graph) const {
	count n = graph.numberOfNodes();
	Permutation permutation;
	std::vector<bool> unvisited(n, true);
	count numUnvisited = n;

	auto remainingArgMinDegree([&](const Graph& graph, const std::vector<bool>& unvisted) {
		index argmin = 0;
		count mindeg = std::numeric_limits<count>::max();

		graph.forNodes([&](node v) {
			if (unvisited.at(v) && graph.degree(v) < mindeg) {
				mindeg = graph.degree(v);
				argmin = v;
			}
		});

		return argmin;
	});


	while (numUnvisited > 0) {
		// for each connected component
		index argmin = remainingArgMinDegree(graph, unvisited); // FIXME: make more efficient
		std::queue<index> q;
		q.push(argmin);
		unvisited.at(argmin) = false;
		--numUnvisited;

		DEBUG("root of RCM: " << argmin);

		while (! q.empty()) {
			// process current front node
			node v = q.front();
			q.pop();
			permutation.push_back(v);

			// enqueue unvisited neighbors in degree-increasing order
			graph.forEdgesOfInDegreeIncreasingOrder(v, [&](node v, node u) {
				if (unvisited.at(u)) {
					q.push(u);
					unvisited.at(u) = false;
					--numUnvisited;
					TRACE("numUnvisited: " << numUnvisited << ", RCM processes " << v);
				}
			});
		}
	}

	// reverse since it is RCM
	std::reverse(permutation.begin(), permutation.end());

	return permutation;
}

} /* namespace NetworKit */
