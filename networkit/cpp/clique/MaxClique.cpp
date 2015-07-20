/*
 * MaxClique.cpp
 *
 *  Created on: 08.12.2014
 *      Author: Henning
 */

#include "MaxClique.h"

namespace NetworKit {

MaxClique::MaxClique(const Graph& G, count lb): G(G), maxi(lb) {

}

void MaxClique::clique(std::set<node>& U, count size) {
	if (U.empty()) {
		if (size > maxi) {
			maxi = size;
		}
	}

	while (U.size() > 0) {
		// pruning 4
		if (size + U.size() <= maxi) {
			return;
		}

		// extract arbitrary element from U
		node x = (* U.begin());
		U.erase(x);

		// pruning 5: compute set of nodes in U that are neighbors of x and have at least degree maxi
		std::set<node> X;
		G.forNeighborsOf(x, [&](node v) {
			if ((G.degree(v) >= maxi) && (U.count(v) > 0)) {
				X.insert(v);
			}
		});

		// recursive call
		clique(X, size + 1);
	}
}

void MaxClique::run() {

	G.forNodes([&](node u) {
		if (G.degree(u) >= maxi) { // pruning 1
			std::set<node> U;
			G.forNeighborsOf(u, [&](node v) {
				if (v > u) { // pruning 2
					if (G.degree(v) >= maxi) { // pruning 3
						U.insert(v);
					}
				}
			});
			clique(U, 1);
		}
	});
}

count MaxClique::getMaxCliqueSize() {
		return maxi;
}

} /* namespace NetworKit */
