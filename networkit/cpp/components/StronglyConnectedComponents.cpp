/*
 * StrongConnectedComponents.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stack>
#include <functional>

#include "StronglyConnectedComponents.h"
#include "../structures/Partition.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

StronglyConnectedComponents::StronglyConnectedComponents(const Graph& G) : G(G) {

}

void StronglyConnectedComponents::run() {
	count z = G.upperNodeIdBound();
	component = Partition(z);

	index nextIndex = 0;
	std::vector<index> nodeIndex(z, none);
	std::vector<index> nodeLowLink(z, none);
	std::stack<node> stx;
	std::vector<bool> onStack(z, false);

	std::function<void(node)> strongConnect = [&](node v) {
		nodeIndex[v] = nextIndex++;
		nodeLowLink[v] = nodeIndex[v];
		stx.push(v);
		onStack[v] = true;

		G.forNeighborsOf(v, [&](node w) {
			if (nodeIndex[w] == none) {
				strongConnect(w);
				nodeLowLink[v] = std::min(nodeLowLink[v], nodeLowLink[w]);
			} else if (onStack[w]) {
				nodeLowLink[v] = std::min(nodeLowLink[v], nodeIndex[w]);
			}
		});

		if (nodeLowLink[v] == nodeIndex[v]) {
			component.toSingleton(v);
			while (true) {
				node w = stx.top();
				stx.pop();
				onStack[w] = false;
				if (w == v) {
					break;
				}
				component[w] = component[v];
			}
		}
	};

	G.forNodes([&](node v) {
		if (nodeIndex[v] == none) {
			strongConnect(v);
		}
	});
}

Partition StronglyConnectedComponents::getPartition() {
	return this->component;
}

count StronglyConnectedComponents::numberOfComponents() {
	return this->component.numberOfSubsets();
}

count StronglyConnectedComponents::componentOfNode(node u) {
	assert (component[u] != none);
	return component[u];
}

}
