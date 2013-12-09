/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiß
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
				if (sl.getShell(w) > i) {
					sl.decreaseShell(w);
				}
			});
		});
	}
	return sl.getShells();

}

} /* namespace NetworKit */
