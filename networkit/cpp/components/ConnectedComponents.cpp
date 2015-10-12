/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include "ConnectedComponents.h"
#include "../auxiliary/Log.h"
#include "../graph/DirOptBFS.h"

namespace NetworKit {

ConnectedComponents::ConnectedComponents(const Graph& G, bool useDirOptBFS) : G(G), hasRun(false), useDirOptBFS(useDirOptBFS) {
	if (G.isDirected()) {
		throw std::runtime_error("Error, connected components of directed graphs cannot be computed, use StronglyConnectedComponents for them.");
	}
}

void ConnectedComponents::run() {
	DEBUG("initializing labels");
	component = Partition(G.upperNodeIdBound(), none);
	numComponents = 0;

	if (!useDirOptBFS) {
		// perform breadth-first searches
		G.forNodes([&](node u) {
			if (component[u] == none) {
				component.setUpperBound(numComponents+1);
				index c = numComponents;
				G.BFSfrom(u, [&](node v) {
					component[v] = c;
				});
				assert (component[u] != none);
				++numComponents;
			}
		});
	} else {
		G.forNodes([&](node u) {
			if (component[u] == none) {
				component.setUpperBound(numComponents+1);
				index c = numComponents;
				DirOptBFS bfs(G, u, false, false);
				std::function<void(node)> callback = [&](node v) {
					component[v] = c;
				};
				bfs.registerCallback(callback);
				bfs.run();
				assert (component[u] != none);
				++numComponents;
			}
		});
	}
	hasRun = true;
}

Partition ConnectedComponents::getPartition() {
	if (!hasRun) throw std::runtime_error("run method has not been called");
	return this->component;
}

std::vector<std::vector<node> > ConnectedComponents::getComponents() {
	if (!hasRun) throw std::runtime_error("run method has not been called");

	// transform partition into vector of unordered_set
	std::vector<std::vector<node> > result(numComponents);

	G.forNodes([&](node u) {
		result[component[u]].push_back(u);
	});

	return result;
}

std::map<index, count> ConnectedComponents::getComponentSizes() {
	if (!hasRun) throw std::runtime_error("run method has not been called");
	return this->component.subsetSizeMap();
}

}
