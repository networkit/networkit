/*
 * RandomSpanningForest.cpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#include "RandomSpanningForest.h"
#include "../graph/Sampling.h"
#include "../components/ConnectedComponents.h"
#include <unordered_set>

namespace NetworKit {

RandomSpanningForest::RandomSpanningForest(const Graph& G):
		SpanningForest(G)
{

}

void RandomSpanningForest::run() {
	// handle disconnected graphs:
	// determine connected components first
	// then start random walk in each component!
	ConnectedComponents cc(G);
	cc.run();
	std::vector<std::vector<node> > comps = cc.getComponents();

	forest = G.copyNodes();
	for (auto comp: comps) {
		std::unordered_set<node> visited;
		const count compSize = comp.size();

		// find and process random root
		index rand = Aux::Random::integer(comp.size() - 1);
		node curr = comp[rand];
		visited.insert(curr);

		// random walk starting from root
		while (visited.size() < compSize) {
			// get random neighbor
			node neigh = G.randomNeighbor(curr);

			// if not seen before, insert tree edge
			if (visited.count(neigh) == 0) {
				forest.addEdge(curr, neigh);
				visited.insert(neigh);
			}

			// move to neighbor
			curr = neigh;
		}
	}
}

} /* namespace NetworKit */
