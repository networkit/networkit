/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include <set>

#include "ConnectedComponents.h"
#include "../structures/Partition.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ConnectedComponents::ConnectedComponents(const Graph& G) : G(G), numComponents(0) {
}

void ConnectedComponents::run() {


	DEBUG("initializing labels");
	component = Partition(G.upperNodeIdBound(), none);
	numComponents = 0;

	// perform breadth-first searches
	G.forNodes([&](node u) {
		if (component[u] == none) {
			component.toSingleton(u);
			index c = component[u];
			assert (component[u] != none);
			G.BFSfrom(u, [&](node v, count dist) {
				component[v] = c;
			});
			++numComponents;
		}
	});

}


Partition ConnectedComponents::getPartition() {
	return this->component;
}


std::map<index, count> ConnectedComponents::getComponentSizes() {
	return this->component.subsetSizeMap();
}

}
