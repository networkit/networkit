/*
 * EffectiveDiameter.cpp
 *
 *  Created on: 16.06.2014
 *      Author: Marc Nemes
 */

#include "EffectiveDiameter.h"
#include "ConnectedComponents.h"

#include <math.h>
#include <iterator>

namespace NetworKit {

count EffectiveDiameter::effectiveDiameter(const Graph& G) {
	// Check whether the graph is connected
	ConnectedComponents cc(G);
	cc.run();
	if (cc.numberOfComponents() > 1) {
		throw std::runtime_error("Graph not connected - diameter is infinite");
	}

	// list of nodes that are connected to all other nodes
	std::set<node> finished;
	// diameter[node][distance][connected_nodes]
	std::vector<std::vector<std::set<node> > > diameter;
	// initialize all nodes
	G.forNodes([&](node v){
		std::set<node> connectedNodes;
		std::vector<std::set<node> > inner;
		// at the beginning, nodes are only connected to themselves
		connectedNodes.insert(v);
		// connect all nodes with themselves at the distance 0
		inner.push_back(connectedNodes);
		diameter.push_back(inner);
	});

	// current diameter
	count d = 1;
	// number of nodes that are connected to all other nodes
	count i = 0;
	// number of nodes that need to be connected with all other nodes
	count threshold = (uint64_t) (ceil(0.9 * G.numberOfNodes()) + 0.5);
	// as long as we need to connect more nodes
	while (i < threshold) {
		G.forNodes([&](node v){
			// only consider nodes that are not already connected to every other node
			if (finished.find(v) == finished.end()) {
				// the node is connected to all nodes from the previous iteration
				std::set<node> connectedNodes = diameter[v][d-1];
				diameter[v].push_back(connectedNodes);
				// and to all previous connected nodes of all neighbors
				G.forNeighborsOf(v, [&](node u) {
					for (auto it : diameter[u][d-1]) {
						// add the current neighbor of u to the neighborhood of v
						diameter[v][d].insert(it);
					}
				});
				// do no longer consider this node once it's connected to all nodes
				if (diameter[v][d].size() == G.numberOfNodes()) {
					finished.insert(v);
					i++;
				}
			}
		});
		d++;
		std::cout << d << std::endl;
	}
	// return the found effective diameter
	return d-1;
}

}
