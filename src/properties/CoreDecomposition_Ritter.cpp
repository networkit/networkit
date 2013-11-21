/*
 * CoreDecomposition_Ritter.cpp
 */

#include "CoreDecomposition_Ritter.h"
#include <list>

namespace NetworKit {

CoreDecomposition_Ritter::CoreDecomposition_Ritter() {

}

CoreDecomposition_Ritter::~CoreDecomposition_Ritter() {

}

std::vector<count> CoreDecomposition_Ritter::run(const Graph& G) {
	const count n = G.numberOfNodes();
	std::vector<count> coreness(n, 0);

	// initializing degree for every node	
	count degree[n];
	count maxDegree = 0;
	G.forNodes([&] (node v) {
		degree[v] = G.degree(v);
		if (degree[v] > maxDegree) {
			maxDegree = degree[v];
		}
	});

	// our main data structure, a list of nodes for each degree
	std::vector< std::list<node> > nodesByDegree(maxDegree + 1);
	// pointer to list entries to move them around fast
	std::list<node>::iterator vptr[n];

	// build initial data structure
	G.forNodes([&](node v) {
		nodesByDegree[degree[v]].push_front(v);
		vptr[v] = nodesByDegree[degree[v]].begin();
	});

	index d_min = 0; // nodesByDegree[d_min] will be our first non empty list
	count i = 0; // current coreness

	for (index j = 0; j < n; j++) {
		// update d_min and i so that nodeByDegree[d_min] is still the first not empty list
		while (nodesByDegree[d_min].empty() && d_min < n) {
			d_min++;
		}
		if (d_min > i) {
			// if we are done the coreness i we can jump to coreness d_min,
			// because if d_min > i there cannot be nodes left with coreness < d_min
			i = d_min;
		}

		// the first node in nodesByDegree[d_min] is the node we are "removing" (and all it edges)
		const node v_remove = nodesByDegree[d_min].front();
		nodesByDegree[d_min].pop_front();
		coreness[v_remove] = i;
		degree[v_remove] = 0; // for remembering that this node is out of the game

		G.forEdgesOf(v_remove, [&] (node v_remove, node w) {
			const index d_old = degree[w];
			const index d_new = d_old - 1;

			if (d_old == 0) { // if the degree is 0, w was already removed from nodesByDegree
				return;
			}

			// update degree and move node to the new list to keep our data structure up to date
			degree[w] = d_new;
			nodesByDegree[d_new].splice(nodesByDegree[d_new].begin(), nodesByDegree[d_old], vptr[w]);

			// may the node w was also in nodesByDegree[d_min] and is now in nodesByDegree[d_min - 1] => adjust d_min
			if (d_new < d_min) {
				d_min = d_new;
			}
		});
	}

	return coreness;
}


} /* namespace NetworKit */
