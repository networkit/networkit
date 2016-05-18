/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß, Christian Staudt
 *  Inplace change on Jun 26, 2015 by Henning Meyerhenke
 */

#include <set>

#include "CoreDecomposition.h"
#include "../auxiliary/PrioQueueForInts.h"
#include <omp.h>

namespace NetworKit {

CoreDecomposition::CoreDecomposition(const Graph& G, bool enforceBucketQueueAlgorithm) :
		Centrality(G, false), maxCore(0), enforceBucketQueueAlgorithm(enforceBucketQueueAlgorithm)
{
	if (G.numberOfSelfLoops()) throw std::runtime_error("Core Decomposition implementation does not support graphs with self-loops. Call Graph.removeSelfLoops() first.");
	canRunInParallel = (! enforceBucketQueueAlgorithm && (G.numberOfNodes() == G.upperNodeIdBound()));
}

void CoreDecomposition::run() {
	if (G.isDirected() || enforceBucketQueueAlgorithm) {
		runWithBucketQueues();
	}
	else {
		runWithParK();
	}
}

void CoreDecomposition::runWithParK() {
	count z = G.upperNodeIdBound();
	scoreData.resize(z); // TODO: move to base class

	count nUnprocessed = G.numberOfNodes();
	std::vector<node> curr; // currently processed nodes
	std::vector<node> next; // nodes to be processed next
	std::vector<char> active(z,0);
	index level = 0; // current level
	count size = 0;  // number of nodes currently processed

	// fill in degrees
	std::vector<count> degrees(z);
	G.parallelForNodes([&](node u) {
		degrees[u] = G.degree(u);
		active[u] = 1;
	});

	// main loop
	while (nUnprocessed > 0) {
		// find nodes with degree == current level
		if (! canRunInParallel || z <= 256) {
			scan(level, degrees, curr);
		}
		else {
			scanParallel(level, degrees, curr, active);
		}

		// process such nodes in curr
		size = curr.size();
		while (size > 0) {
			nUnprocessed -= size;
			if (! canRunInParallel || size <= 256) {
				processSublevel(level, degrees, curr, next);
			}
			else {
				processSublevelParallel(level, degrees, curr, next, active);
			}

			std::swap(curr, next);
			size = curr.size();
			next.clear();
		}
		++level;
	}

	maxCore = level-1;
	hasRun = true;
}

void NetworKit::CoreDecomposition::scan(index level, const std::vector<count>& degrees,
		std::vector<node>& curr)
{
	G.forNodes([&](node u) {
		if (degrees[u] == level) {
			curr.push_back(u);
		}
	});
}

void NetworKit::CoreDecomposition::scanParallel(index level, const std::vector<count>& degrees,
		std::vector<node>& curr, std::vector<char>& active)
{
	const count z = G.upperNodeIdBound();
	std::vector<std::vector<node>> next(omp_get_max_threads());
	curr.clear();

#pragma omp parallel for schedule(guided)
	for (index u = 0; u < z; ++u) {
		if (active[u] && degrees[u] == level) {
			auto tid = omp_get_thread_num();
			next[tid].push_back(u);
		}
	}
	for (auto& n : next) {
		curr.insert(curr.end(),n.begin(),n.end());
	}
}

void NetworKit::CoreDecomposition::processSublevel(index level,
		std::vector<count>& degrees, const std::vector<node>& curr,
		std::vector<node>& next)
{
	// check for each neighbor of vertices in curr if their updated degree reaches level;
	// if so, process them next
	for (auto u: curr) {
		scoreData[u] = level;
		G.forNeighborsOf(u, [&](node v) {
			if (degrees[v] > level) {
				degrees[v]--;
				if (degrees[v] == level) {
					next.push_back(v);
				}
			}
		});
	}
}

void NetworKit::CoreDecomposition::processSublevelParallel(index level,
		std::vector<count>& degrees, const std::vector<node>& curr,
		std::vector<node>& next, std::vector<char>& active)
{
	// check for each neighbor of vertices in curr if their updated degree reaches level;
	// if so, process them next

	const count size = curr.size();
	std::vector<std::vector<node>> localNext(omp_get_max_threads());

#pragma omp parallel for schedule(guided)
	for (index i = 0; i < size; ++i) {
		node u = curr[i];
		active[u] = 0;
		scoreData[u] = level;
		G.forNeighborsOf(u, [&](node v) {
			if (degrees[v] > level) {
				index tmp;
#pragma omp atomic capture
				tmp = --degrees[v];

				// ensure that neighbor is inserted exactly once if necessary
				if (tmp == level) { // error in external publication on ParK fixed here
					auto tid = omp_get_thread_num();
					localNext[tid].push_back(v);
				}
			}
		});
	}
	for (auto& n : localNext) {
		next.insert(next.end(), n.begin(), n.end());
	}
}

void CoreDecomposition::runWithBucketQueues() {
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


Cover CoreDecomposition::getCover() const {
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

Partition CoreDecomposition::getPartition() const {
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

double CoreDecomposition::maximum() {
	return G.numberOfNodes() - 1;
}

} /* namespace NetworKit */

