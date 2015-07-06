/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß, Christian Staudt
 *  Inplace change on Jun 26, 2015 by Henning Meyerhenke
 */

#include <set>

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition(const Graph& G) : Centrality(G, false), maxCore(0) {
	if (G.numberOfSelfLoops()) throw std::runtime_error("Core Decomposition implementation does not supprt graphs with self-loops. Call Graph.removeSelfLoops() first.");
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

	std::vector<count> degree(z);       // tracks degree during algo
	std::vector<bool> alive(z, true);   // tracks if node is deleted (alive) or not
	count numAlive = G.numberOfNodes(); // tracks number of alive nodes

	/* Insert nodes into their initial buckets. */
	if (!G.isDirected()) {
		G.forNodes([&](node v) {
			count deg = G.degree(v);
			degree[v] = deg;
			buckets[deg].push_front(v);
			core = std::min(core, deg);
			nodePtr[v] = buckets[deg].begin();
		});
	} else {
		G.forNodes([&](node v) {
			count deg = G.degreeIn(v) + G.degreeOut(v); // TODO: Document this behavior for directed graph
			degree[v] = deg;
			buckets[deg].push_front(v);
			core = std::min(core, deg);
			nodePtr[v] = buckets[deg].begin();
		});

	}

	/* Main loop: Successively "remove" nodes by setting them not alive after processing them. */
	while (numAlive > 0) {
		Bucket& cur_bucket = buckets[core];

		/* Remove nodes with remaining degree <= core. */
		while (!cur_bucket.empty()) {
			/* scoreData for node u is current core value. */
			node u = cur_bucket.front();
			cur_bucket.pop_front();
			scoreData[u] = core;

			/* Remove u and its incident edges. */
			/* graph is undirected */
			if (!G.isDirected()) {
				G.forNeighborsOf(u, [&](node v) {
					if (alive[v]) {
						count deg = degree[v];
						degree[v]--;

						/* Shift node v into new bucket.
					   Optimisation: Need not move to buckets < core. */
						if (deg > core) {
							buckets[deg].erase(nodePtr[v]);
							buckets[deg - 1].push_front(v);
							nodePtr[v] = buckets[deg - 1].begin();
						}
					}
				});
			} else {
			/* graph is directed */
				G.forNeighborsOf(u, [&](node v) {
					if (alive[v]) {
						count deg = degree[v];
						degree[v]--;

						/* Shift node v into new bucket.
					   Optimisation: Need not move to buckets < core. */
						if (deg > core) {
							buckets[deg].erase(nodePtr[v]);
							buckets[deg - 1].push_front(v);
							nodePtr[v] = buckets[deg - 1].begin();
						}
					}
				});
				G.forInNeighborsOf(u, [&](node u, node v) {
					if (alive[v]) {
						count deg = degree[v];
						degree[v]--;

						/* Shift node v into new bucket.
					   Optimisation: Need not move to buckets < core. */
						if (deg > core) {
							buckets[deg].erase(nodePtr[v]);
							buckets[deg - 1].push_front(v);
							nodePtr[v] = buckets[deg - 1].begin();
						}
					}
				});
			}

			// "delete" current node
			alive[u] = false;
			--numAlive;
		}
		core++;
	}

	maxCore = core - 1;

	hasRun = true;
}


Cover CoreDecomposition::cores() const {
	if (! hasRun) throw std::runtime_error("call run method first");
	// initialize Cover
	index z = G.upperNodeIdBound();
	Cover coverData;
	if (coverData.numberOfElements() != z) {
		coverData = Cover(z);
		coverData.setUpperBound(z);
	}
	// enter values from scoreData into coverData
	G.forNodes([&](node u) {
		index k = 0;
		while (scoreData[u] >= k) {
			coverData.addToSubset((index) k, (index) u);
			++k;
		}
	});
	return coverData;
}

Partition CoreDecomposition::shells() const {
	if (! hasRun) throw std::runtime_error("call run method first");
	// initialize Partition
	index z = G.upperNodeIdBound();
	Partition shellData;
	if (shellData.numberOfElements() != z) {
		shellData = Partition(z);
		shellData.allToSingletons();
	}
	// enter values from scoreData into shellData
	G.forNodes([&](node u){
		shellData.moveToSubset((index) scoreData[u], (index) u);
	});
	return shellData;
}

index CoreDecomposition::maxCoreNumber() const {
	if (! hasRun) throw std::runtime_error("call run method first");
	return maxCore;
}

} /* namespace NetworKit */
