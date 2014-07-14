/*
 * RegionGrowingOverlapper.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RegionGrowingOverlapper.h"
#include "../auxiliary/Log.h"

namespace NetworKit {



Partition RegionGrowingOverlapper::run(const Graph& G,
		const std::vector<Partition>& clusterings) {
	// TODO: test

	uint64_t n = G.numberOfNodes();
	Partition core(n); // "core groups" resulting from overlap of all clustering
	core.allToSingletons(); // assign all nodes to singletons

	std::vector<int> visited(n, 0); // node -> has been visited (1) or not (0). not <bool> because of thread-safety

	node r = -1; // start node for BFS
	//bool allVisited = false; // have all nodes been visited by BFS?

	std::set<index> unvisited;
	for (count i = 0; i < n; ++i) {
		unvisited.insert(i);
	}

	while (unvisited.size() > 0) { // !allVisited) {
		// select start node for BFS

//		for (node v = 0; v < n; ++v) {
//			if (visited[v] == 0) {
//				r = v;
//				break;
//			}
//			allVisited = true;
//		}

		// root: first unvisited element
		r = *unvisited.begin();
		unvisited.erase(r);

		DEBUG("starting BFS from node: " , r);

		// TODO: Graph.breadthFirstNodesFrom() was removed, use BFSfrom instead - marvin
		// <<<<< old code
		// // start BFS
		// G.breadthFirstNodesFrom(r, visited, [&](node u) {
		// 	visited[u] = 1; // has been visited
		// 		unvisited.erase(u);
		// 		// check for all incident edges if u and v belong in the same core cluster
		// 		G.forEdgesOf(u, [&](node u, node v) {
		// 					bool together = true;
		// 					for (std::vector<Partition>::iterator iter = clusterings.begin(); iter != clusterings.end(); ++iter ) {
		// 						together = together && (iter->subsetOf(u) == iter->subsetOf(v));
		// 					}
		// 					if (together) {
		// 						core.moveToSubset(core.subsetOf(u), v);
		// 					}
		// 				});
		// 	});
		// >>>>>

		// <<<<< new code
		G.BFSfrom(r, [&](node u) {
			unvisited.erase(u);
			// check for all incident edges if u and v belong in the same core cluster
			G.forEdgesOf(u, [&](node u, node v) {
				bool together = true;
				for (auto iter = clusterings.begin(); iter != clusterings.end(); ++iter ) {
					together = together && (iter->subsetOf(u) == iter->subsetOf(v));
				}
				if (together) {
					core.moveToSubset(core.subsetOf(u), v);
				}
			});
		});
		// >>>>>
	}

	return core;
}

} /* namespace NetworKit */
