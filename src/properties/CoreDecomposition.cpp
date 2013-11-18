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
	std::vector<count> coreness;

	Aux::ShellList sl(&G);

	G.forNodes([&](node v) {
		sl.insert(v);
	});

	for (int i = 0; i < sl.size(); i++) {
		std::list<node>::iterator it = sl.getShelliterator(i);

		sl.forEachNodeInShell(i, [&](node n) {
			G.forNeighborsOf(n, [&](node m) {
				if (sl.getCurrentShell(m) > i) {
					sl.decreaseDegree(m);
				}
			});
		});
	}

	this->coreness = sl.getCoreness();
	return this->coreness;

}

void
CoreDecomposition::write(std::string filename) {
	std::ofstream out;
  	out.open(filename);

  	for (auto it = this->coreness.begin(); it != this->coreness.end(); ++it) {
  		out << *it << std::endl;
  	}

  	out.close();
 }

} /* namespace NetworKit */
