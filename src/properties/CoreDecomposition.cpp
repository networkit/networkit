/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition() {

}

CoreDecomposition::~CoreDecomposition() {

}

std::vector<count> CoreDecomposition::run(const Graph& G) {
	Aux::ShellList sl(&G);

	for (count i = 0; i < sl.size(); i++) {
		sl.forEachNodeInShell(i, [&](node v) {
			G.forNeighborsOf(v, [&](node w) {
				if (sl.getCurrentShell(w) > i) {
					sl.decreaseShell(w);
				}
			});
		});
	}
	return sl.getCoreness();

}

} /* namespace NetworKit */
