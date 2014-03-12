/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include "ConnectedComponents.h"
#include "../structures/Partition.h"
#include "../coarsening/PartitionCoarsening.h"
#include <set>

namespace NetworKit {

ConnectedComponents::ConnectedComponents() {

}

ConnectedComponents::~ConnectedComponents() {

}

void ConnectedComponents::run(const Graph& G) {
	// calculate connected components by label propagation
	count z = G.numberOfNodes();
	this->component = std::vector<node>(z, none);

	DEBUG("initializing labels");
	Partition partition(G.upperNodeIdBound());
	partition.allToSingletons();

	DEBUG("initializing active nodes");
	std::vector<bool> activeNodes(z); // record if node must be processed
	activeNodes.assign(z, true);
	G.forNodes([&](node u) { // NOTE: not in parallel due to implementation of bit vector
		if (G.degree(u) == 0) {
			activeNodes[u] = false;
		}
	});

	DEBUG("main loop");
//	count numActive = 0; // for debugging purposes only
	count numIterations = 0;
	bool change = false;
	do {
//		TRACE("label propagation iteration");
		change = false;
//		numActive = 0;
		G.forNodes([&](node u) {
			if (activeNodes[u]) {
//				++numActive;
				std::vector<index> neighborLabels;
				G.forNeighborsOf(u, [&](node v) {
					// neighborLabels.push_back(component[v]);
					neighborLabels.push_back(partition[v]);
				});
				// get smallest
				index smallest = *std::min_element(neighborLabels.begin(), neighborLabels.end());

				if (partition[u] != smallest) {
					partition.moveToSubset(smallest, u);
					change = true;
					G.forNeighborsOf(u, [&](node v) {
						activeNodes[v] = true;
					});
				} else {
					activeNodes[u] = false; // current node becomes inactive
				}
			}
		});
//		TRACE("num active: ", numActive);
		++numIterations;
		if ((numIterations % 8) == 0) { // TODO: externalize constant
			// coarsen and make recursive call
			PartitionCoarsening con;
			std::pair<Graph, std::vector<node> > coarse = con.run(G, partition);
			ConnectedComponents cc;
			cc.run(coarse.first);

			// apply to current graph
			G.forNodes([&](node u) {
				partition[u] = cc.componentOfNode(coarse.second[u]);
			});
		}
	} while (change);

	G.parallelForNodes([&](node u) {
		component[u] = partition[u];
	});
}


std::vector<index> ConnectedComponents::getComponentData() {
	return this->component;
}

std::map<index, count> ConnectedComponents::getComponentSizes() {
	std::map<index, count> componentSize;
	for (node u = 0; u < component.size(); ++u) {
		if (component[u] != none) {
			componentSize[component[u]] += 1;
		}
	}
	return componentSize;
}


std::vector<node> ConnectedComponents::getComponent(index label) {
	std::vector<node> nodesOfComponent;
	for (node u = 0; u < component.size(); ++u) {
		if (this->component[u] == label) {
			nodesOfComponent.push_back(u);
		}
	}
	return nodesOfComponent;
}

count ConnectedComponents::numberOfComponents() {
	count nc = 0;
	std::set<index> componentLabels;
	for (index i = 0; i < component.size(); ++i) {
		index c = component[i];
		if ((componentLabels.find(c) == componentLabels.end()) && (c != none)) { // label not encountered yet
			componentLabels.insert(c);
			nc++;
		}
	}
	return nc;
}

count ConnectedComponents::componentOfNode(node u) {
	assert (component[u] != none);
	return component[u];
}

}

