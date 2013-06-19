/*
 * PseudoDynamic.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "PseudoDynamic.h"

namespace NetworKit {

PseudoDynamic::PseudoDynamic(const Graph& Gstatic) : Gstatic(Gstatic), u(0) {
	// TODO Auto-generated constructor stub

}

PseudoDynamic::~PseudoDynamic() {
	// TODO Auto-generated destructor stub
}

void PseudoDynamic::initializeGraph() {
}

void PseudoDynamic::generate() {
	if (this->G->numberOfNodes() == Gstatic.numberOfNodes()) {
		throw std::logic_error("all nodes from the static graph have been generated");
	}
	u += 1;
	if (Gstatic.hasNode(u)) {
		node uNew = this->Gproxy->addNode();
		assert (u == uNew);

		Gstatic.forNeighborsOf(u, [&](node v){
			if (v <= u) {
				this->Gproxy->addEdge(u, v);
			}
		});
	}
}

} /* namespace NetworKit */
