/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß, Christian Staudt
 */

#include <set>

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition(const Graph& G) : Centrality(G, false), maxCore(0) {

}

void CoreDecomposition::run() {
	/* Main data structure: buckets of nodes indexed by their remaining degree. */
	typedef std::list<node> Bucket;
	index z = G.upperNodeIdBound();
	std::vector<Bucket> buckets(z);
	std::vector<Bucket::iterator> nodePtr(z);

	/* Current core and and computed scoreData values. */
	index core = std::numeric_limits<index>::max();
	scoreData.clear();
	scoreData.resize(z);

	/* Insert nodes into their initial buckets. */
	if (!G.isDirected()) {
		G.forNodes([&](node v) {
			count deg = G.degree(v);
			buckets[deg].push_front(v);
			core = std::min(core, deg);
			nodePtr[v] = buckets[deg].begin();
		});
	} else {
		G.forNodes([&](node v) {
			count deg = G.degreeIn(v) + G.degreeOut(v);
			buckets[deg].push_front(v);
			core = std::min(core, deg);
			nodePtr[v] = buckets[deg].begin();
		});

	}

	/* Main loop: Successively remove nodes in copy G2 of G. */
	Graph G2 = G;
	while (!G2.isEmpty()) {
		Bucket& cur_bucket = buckets[core];

		/* Remove nodes with remaining degree <= core. */
		while (!cur_bucket.empty()) {
			/* scoreData for node u is current core value. */
			node u = cur_bucket.front();
			cur_bucket.pop_front();
			scoreData[u] = core;

			/* Remove u and its incident edges. */
			/* graph is undirected */
			if (!G2.isDirected()) {
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
			} else {
			/* graph is directed */
				G2.forNeighborsOf(u, [&](node v) {
					count deg = G2.degreeIn(v) + G2.degreeOut(v);
					G2.removeEdge(u, v);

					/* Shift node v into new bucket.
					   Optimisation: Need not move to buckets < core. */
					if (deg > core) {
						buckets[deg].erase(nodePtr[v]);
						buckets[deg - 1].push_front(v);
						nodePtr[v] = buckets[deg - 1].begin();
					}
				});
				G2.forInNeighborsOf(u, [&](node u, node v) {
					count deg = G2.degreeOut(v) + G2.degreeIn(v);
					G2.removeEdge(v, u);

					/* Shift node v into new bucket.
					   Optimisation: Need not move to buckets < core. */
					if (deg > core) {
						buckets[deg].erase(nodePtr[v]);
						buckets[deg - 1].push_front(v);
						nodePtr[v] = buckets[deg - 1].begin();
					}
				});
			}
			G2.removeNode(u);
		}
		core++;
	}

	maxCore = core - 1;
	hasRun = true;
}


std::vector<std::set<node> > CoreDecomposition::cores() const {
	if (! hasRun) throw std::runtime_error("call run method first");

	std::vector<std::set<node> > cores(maxCore + 1);
	for (index k = 0; k <= maxCore; k++) {
		G.forNodes([&](node u){
			if (scoreData[u] >= k) {
				cores.at(k).insert(u);
			}
		});
	}
	return cores;
}

std::vector<std::set<node> > CoreDecomposition::shells() const {
	if (! hasRun) throw std::runtime_error("call run method first");

	std::vector<std::set<node> > shells(maxCore + 1);
	for (index k = 0; k <= maxCore; k++) {
		G.forNodes([&](node u){
			if (scoreData[u] == k) {
				shells.at(k).insert(u);
			}
		});
	}
	return shells;
}

index CoreDecomposition::maxCoreNumber() const {
	if (! hasRun) throw std::runtime_error("call run method first");
	return maxCore;
}

} /* namespace NetworKit */
