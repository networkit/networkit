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
	// do nothing - start with an empty graph
}

void PseudoDynamic::generate() {
	TRACE("calling generate");
	// if there are no more nodes left from the input graph, stop
	if (this->G->numberOfNodes() == Gstatic.numberOfNodes()) {
		throw std::logic_error("all nodes from the static graph have been generated");
	}


	if (Gstatic.hasNode(u)) {
		TRACE("static graph has node " << u);
		// create node in the dynamic graph
		node uNew = this->Gproxy->addNode();
		TRACE("creating new node " << uNew);

		// add edges to already existing neighbors
		Gstatic.forNeighborsOf(u, [&](node v){
			if (v <= u) { // assumption: all nodes below the current one exist in the dynamic graph
				this->Gproxy->addEdge(u, v);
			}
		});
	} else {
		TRACE("node does not exist in the static graph: " << u);
	}

	// next node
	u += 1;
	// trigger time step
	this->Gproxy->timeStep();
}

} /* namespace NetworKit */
