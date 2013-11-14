/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning, Barth, Weiß
 */

#include "../auxiliary/ShellList.h"
#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition() {
}

CoreDecomposition::~CoreDecomposition() {
}

std::vector<count> CoreDecomposition::run(const Graph& G) {
	Aux::ShellList sl(&G);

  /* Main Loop
   * Iterates over all shells beginning with the first. The zero shell can safely be ignored!
   */
  for (count i = 1; i < sl.size(); i++) {
    sl.forEachNodeInShell(i, [&](node v) {
      /* Shell Nodes Loop
       * Within each shell, all neighbors of the current node fall down one shell if they are in a higher one.
       */
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
