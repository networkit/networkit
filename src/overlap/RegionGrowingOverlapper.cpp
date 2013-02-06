/*
 * RegionGrowingOverlapper.cpp
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#include "RegionGrowingOverlapper.h"

namespace EnsembleClustering {

RegionGrowingOverlapper::RegionGrowingOverlapper() {
	// TODO Auto-generated constructor stub

}

RegionGrowingOverlapper::~RegionGrowingOverlapper() {
	// TODO Auto-generated destructor stub
}

Clustering RegionGrowingOverlapper::run(Graph& G, std::vector<Clustering>& clusterings) {
	// TODO: test

	int64_t n = G.numberOfNodes();
	Clustering core(n); // "core groups" resulting from overlap of all clustering
	core.allToSingletons(); // assign all nodes to singletons

	NodeMap<int> visited(n, 0); // node -> has been visited (1) or not (0). not <bool> because of thread-safety
	node r;	// start node for BFS
	bool allVisited = false;	// have all nodes been visited by BFS?

	while (! allVisited) {
		// select start node for BFS

		for (node v = 0; v < n; ++v) {
			if (visited[v] == 0) {
				r = v;
				break;
			}
			allVisited = true;
		}

		// start BFS
		G.breadthFirstNodesFrom(r, [&](node u) {
			visited[u] = 1; // has been visited

			// check for all incident edges if u and v belong in the same core cluster
			G.forEdgesOf(u, [&](node u, node v) {
				bool together = true;
				// TODO: is this cache-efficient? if not: use combined clustering array with k entries per node
				for (Clustering zeta : clusterings) {
					together = together && (zeta.clusterOf(u) == zeta.clusterOf(v));
				}
				if (together) {
					core.moveToCluster(core.clusterOf(u), v);
				}
			});

		});
	}

	return core;
}

} /* namespace EnsembleClustering */
