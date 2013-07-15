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


std::map<index, index> RcmMapper::run(Graph& guest, Graph& host) {
	std::map<index, index> mapping;


	// TODO: permute guest graph according to RCM

	// TODO: permute host graph according to RCM

	// TODO: fill in mapping according to permutations


	return mapping;
}

std::vector<index> RcmMapper::permute(const Graph& graph) const {
	std::vector<index> permutation;
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
			graph.forNeighborsOf(v, [&](node u) { // FIXME: InDegreeIncreasingOrder(v, [&](node u) {
				if (unvisited.count(u) > 0) {
					q.push(u);
				}
			});
		}
	}

	// TODO: perform RCM

	return permutation;
}

} /* namespace NetworKit */
