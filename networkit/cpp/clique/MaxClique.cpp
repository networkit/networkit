/*
 * MaxClique.cpp
 *
 *  Created on: 08.12.2014
 *      Author: Henning
 */

#include "MaxClique.h"
#include "../auxiliary/SignalHandling.h"


namespace NetworKit {

MaxClique::MaxClique(const Graph& G, count lb): G(G), maxi(lb) {

}

void MaxClique::clique(std::unordered_set<node>& U, std::unordered_set<node>& currClique, count size) {
	Aux::SignalHandler handler;

	if (U.empty()) {
		if (size > maxi) {
			maxi = size;
			bestClique = currClique;
			INFO("new best clique, size: ", maxi);
//			assert(size == bestClique.size());
		}
		return;
	}

	while (U.size() > 0) {
		if (! handler.isRunning()) {
			return;
		}

		// pruning 4
		if (size + U.size() <= maxi) {
			return;
		}

		// extract arbitrary element from U
		node x = (* U.begin());
		U.erase(x);

		// pruning 5: compute unordered_set of nodes in U that are neighbors of x and have at least degree maxi
		std::unordered_set<node> X;
		G.forNeighborsOf(x, [&](node v) {
			if ((G.degree(v) >= maxi) && (U.count(v) > 0)) {
				X.insert(v);
			}
		});

		// recursive call
		std::unordered_set<node> extendedClique = currClique;
		extendedClique.insert(x);
//		assert(extendedClique.size() == size + 1);
		clique(X, extendedClique, size + 1);
	}
}

void MaxClique::run() {
	Aux::SignalHandler handler;
	std::unordered_set<node> currClique;

	G.forNodes([&](node u) {
		if (! handler.isRunning()) {
			return;
		}

		if (G.degree(u) >= maxi) { // pruning 1
			std::unordered_set<node> U;
			G.forNeighborsOf(u, [&](node v) {
				if (v > u) { // pruning 2
					if (G.degree(v) >= maxi) { // pruning 3
						U.insert(v);
					}
				}
			});
			currClique.clear();
			currClique.insert(u);
//			assert(currClique.size() == 1);
//			if (u % 4 == 0)
//				INFO("initial recursion for ", u);
			clique(U, currClique, 1);
		}
	});
}

count MaxClique::getMaxCliqueSize() {
	return maxi;
}

std::unordered_set<node> MaxClique::getMaxClique() const {
	return bestClique;
}

} /* namespace NetworKit */
