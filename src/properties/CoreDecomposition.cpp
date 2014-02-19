/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß
 */

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition() {
}

CoreDecomposition::~CoreDecomposition() {
}

std::vector<count> CoreDecomposition::run(const Graph& G) {
	/* Main data structure: buckets of nodes indexed by their remaining degree. */
	typedef std::list<node> Bucket;
	count nnodes = G.numberOfNodes();
	std::vector<Bucket> buckets(nnodes);
	std::vector<Bucket::iterator> nodePtr(nnodes);

	/* Current core and and computed coreness values. */
	count core = std::numeric_limits<count>::max();
	std::vector<count> coreness(nnodes);

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

	return coreness;
}

} /* namespace NetworKit */
