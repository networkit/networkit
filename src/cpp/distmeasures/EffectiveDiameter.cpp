#include "EffectiveDiameter.h"

namespace NetworKit {
	int EffectiveDiameter::effectiveDiameter(const Graph& G) {
		// list of nodes that are connected to all other nodes
		std::set<int> finished;
		// diameter[node][distance][connected_nodes]
		std::vector<std::vector<std::set<int> > > diameter (G.numberOfNodes()); //
		// initialize all nodes
		G.forNodes([&](node v){
			std::set<int> list;
			// at the beginning, nodes are only connected to themselves
			list.insert(v);
			// connect all nodes with themselves at the distance 0
			diameter[v][0] = list;
		});

		// current diameter
		int d = 1;
		// number of nodes that are connected to all other nodes
		uint64_t i = 0;
		// number of nodes that need to be connected with all other nodes
		uint64_t threshold = (uint64_t) (ceil(0.9 * G.numberOfNodes()) + 0.5);
		// as long as we need to connect more nodes
		while (i < threshold) {
			G.forNodes([&](node v){
				// only consider nodes that are not already connected to every other node
				if (finished.find(v) == finished.end()) {
					// the node is connected to all nodes from the previous iteration
					diameter[v][d] = diameter[v][d-1];
					// and to all previous connected nodes of all neighbors
					G.forNeighborsOf(v, [&](node u) {
						//diameter[v][d] = diameter[u][d-1];
						std::set<int>::iterator it;
						for (it = diameter[u][d-1].begin(); it != diameter[u][d-1].end(); ++it)	{
						    // add the current neighbor of u to the neighborhood of v
						    diameter[v][d].insert(*it);
						}
					});
					// do no longer consider this node once it's connected to all nodes
					if (diameter[v][d].size() == G.numberOfNodes()) {
						finished.insert(v);
						i++;
					}
				}
			});
		}
		// return the found effective diameter
		return d;
	}
}
