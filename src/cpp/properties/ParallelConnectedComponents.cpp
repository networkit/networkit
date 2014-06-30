/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include <set>

#include "ParallelConnectedComponents.h"
#include "../structures/Partition.h"
#include "../coarsening/PartitionCoarsening.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ParallelConnectedComponents::ParallelConnectedComponents(const Graph& G, bool coarsening) : G(G), coarsening(coarsening) {

}


void ParallelConnectedComponents::run() {
	if (G.isDirected()) {
		throw std::runtime_error("algorithm does not accept directed graphs");
	}

	// calculate connected components by label propagation
	count z = G.numberOfNodes();

	DEBUG("initializing labels");
	component = Partition(G.upperNodeIdBound());
	component.allToSingletons();

	DEBUG("initializing active nodes");
	const char INACTIVE = 0;
	const char ACTIVE = 1;
	std::vector<char> activeNodes(z); // record if node must be processed
	std::vector<char> nextActiveNodes(z, ACTIVE); // for next iteration
	nextActiveNodes.assign(z, ACTIVE);
	G.forNodes([&](node u) { // NOTE: not in parallel due to implementation of bit vector
		if (G.degree(u) == 0) {
			nextActiveNodes[u] = INACTIVE;
		}
	});

	DEBUG("main loop");
//	count numActive = 0; // for debugging purposes only
	count numIterations = 0;
	bool change = false;
	do {
//		TRACE("label propagation iteration");
		activeNodes = nextActiveNodes;
		nextActiveNodes.assign(z, INACTIVE);
		change = false;
//		numActive = 0;
		G.parallelForNodes([&](node u) {
			if (activeNodes[u] == ACTIVE) {
//				++numActive;
				// get smallest
				index smallest = component[u];
				G.forNeighborsOf(u, [&](node v) {
					smallest = std::min(smallest, component[v]);
				});

				if (component[u] != smallest) {
					component.moveToSubset(smallest, u);
					change = true;
					G.forNeighborsOf(u, [&](node v) {
						nextActiveNodes[v] = ACTIVE;
					});
				}
//				else {
//					nextActiveNodes[u] = INACTIVE; // current node becomes inactive
//				}
			}
		});
//		TRACE("num active: ", numActive);
		++numIterations;
		if (coarsening && (numIterations % 8) == 0) { // TODO: externalize constant
			// coarsen and make recursive call
			PartitionCoarsening con;
			std::pair<Graph, std::vector<node> > coarse = con.run(G, component);
			ParallelConnectedComponents cc(coarse.first);
			cc.run();

			// apply to current graph
			G.forNodes([&](node u) {
				component[u] = cc.componentOfNode(coarse.second[u]);
			});
		}
	} while (change);
}


void ParallelConnectedComponents::runSequential() {
	if (G.isDirected()) {
		throw std::runtime_error("algorithm does not accept directed graphs");
	}
	// calculate connected components by label propagation
	count z = G.numberOfNodes();
	DEBUG("initializing labels");
	component = Partition(G.upperNodeIdBound());
	component.allToSingletons();
	DEBUG("initializing active nodes");
	std::vector<bool> activeNodes(z); // record if node must be processed
	activeNodes.assign(z, true);
	G.forNodes([&](node u) { // NOTE: not in parallel due to implementation of bit vector
		if (G.degree(u) == 0) {
			activeNodes[u] = false;
		}
	});
	DEBUG("main loop");
	// count numActive = 0; // for debugging purposes only
	count numIterations = 0;
	bool change = false;
	do {
		// TRACE("label propagation iteration");
		change = false;
		// numActive = 0;
		G.forNodes([&](node u) {
			if (activeNodes[u]) {
				// ++numActive;
				std::vector<index> neighborLabels;
				G.forNeighborsOf(u, [&](node v) {
					// neighborLabels.push_back(component[v]);
					neighborLabels.push_back(component[v]);
				});
				// get smallest
				index smallest = *std::min_element(neighborLabels.begin(), neighborLabels.end());
				if (component[u] != smallest) {
					component.moveToSubset(smallest, u);
					change = true;
					G.forNeighborsOf(u, [&](node v) {
						activeNodes[v] = true;
					});
				} else {
					activeNodes[u] = false; // current node becomes inactive
				}
			}
		});
		// TRACE("num active: ", numActive);
		++numIterations;
		if ((numIterations % 8) == 0) { // TODO: externalize constant
			// coarsen and make recursive call
			PartitionCoarsening con;
			std::pair<Graph, std::vector<node> > coarse = con.run(G, component);
			ParallelConnectedComponents cc(coarse.first);
			cc.run();
			// apply to current graph
			G.forNodes([&](node u) {
				component[u] = cc.componentOfNode(coarse.second[u]);
			});
		}
	} while (change);
}


Partition ParallelConnectedComponents::getPartition() {
	return this->component;
}

count ParallelConnectedComponents::numberOfComponents() {
	return this->component.numberOfSubsets();
}

count ParallelConnectedComponents::componentOfNode(node u) {
	assert (component[u] != none);
	return component[u];
}

}
