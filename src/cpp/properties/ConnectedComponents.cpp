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

ConnectedComponents::ConnectedComponents(const Graph& G) : G(G) {

}


void ConnectedComponents::run() {


	DEBUG("initializing labels");
	component = Partition(G.upperNodeIdBound(), none);

	// perform breadth-first searches
	G.forNodes([&](node u) {
		TRACE("scanning node ", u);
		if (component[u] == none) {
			TRACE(u, "is not assigned to component");
			component.toSingleton(u);
			index c = component[u];
			TRACE("new component of ", u, ": ", c);
			assert (component[u] != none);
			G.BFSfrom(u, [&](node v) {
				component[v] = c;
				TRACE("assigned ", v, " to ", c);
			});
		}
	});

}


Partition ConnectedComponents::getPartition() {
	return this->component;
}

count ConnectedComponents::numberOfComponents() {
	return this->component.numberOfSubsets();
}

count ConnectedComponents::componentOfNode(node u) {
	assert (component[u] != none);
	return component[u];
}

}
