/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß, Christian Staudt
 */

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition(const Graph& G) : G(G) {

}

void CoreDecomposition::run() {
	/* Main data structure: buckets of nodes indexed by their remaining degree. */
	typedef std::list<node> Bucket;
	index z = G.upperNodeIdBound();
	std::vector<Bucket> buckets(z);
	std::vector<Bucket::iterator> nodePtr(z);

	/* Current core and and computed coreness values. */
	index core = std::numeric_limits<index>::max();
	coreness.clear();
	coreness.resize(z);

	/* Insert nodes into their initial buckets. */
	G.forNodes([&](node v) {
		count deg = G.degree(v);
		buckets[deg].push_front(v);
		core = std::min(core, deg);
		nodePtr[v] = buckets[deg].begin();
	});

	/* Main loop: Successively remove nodes in copy G2 of G. */
	Graph G2 = G;
	while (!G2.isEmpty()) {
		Bucket& cur_bucket = buckets[core];

		/* Remove nodes with remaining degree <= core. */
		while (!cur_bucket.empty()) {
			/* Coreness for node u is current core value. */
			node u = cur_bucket.front();
			cur_bucket.pop_front();
			coreness[u] = core;

			/* Remove u and its incident edges. */
			G2.forNeighborsOf(u, [&](node v) {
				count deg = G2.degree(v);
				G2.removeEdge(u, v);

				/* Shift node v into new bucket.
				   Optimisation: Need not move to buckets < core. */
				if (deg > core) {
					buckets[deg].erase(nodePtr[v]);
					buckets[deg - 1].push_front(v);
					nodePtr[v] = buckets[deg - 1].begin();
				}
			});
			G2.removeNode(u);
		}
		core++;
	}

	maxCoreNumber = core - 1;
}

std::vector<index> CoreDecomposition::coreNumbers() const {
	return coreness;
}

index CoreDecomposition::coreNumber(node v) const {
	return coreness.at(v);
}


std::vector<std::set<node> > CoreDecomposition::cores() const {
	std::vector<std::set<node> > cores(maxCoreNumber);
	for (index k = 0; k <= maxCoreNumber; k++) {
		G.forNodes([&](node u){
			if (coreness[u] >= k) {
				cores.at(k).insert(u);
			}
		});
	}
	return cores;
}

std::vector<std::set<node> > CoreDecomposition::shells() const {
	std::vector<std::set<node> > shells;
	for (index k = 0; k <= maxCoreNumber; k++) {
		G.forNodes([&](node u){
			if (coreness[u] == k) {
				shells.at(k).insert(u);
			}
		});
	}
	return shells;
}

} /* namespace NetworKit */
