/*
 * PseudoDynamic.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "PseudoDynamic.h"

namespace NetworKit {

PseudoDynamic::PseudoDynamic(Graph& Gstatic) : Gstatic(Gstatic), idmap(Gstatic.upperNodeIdBound(), none), current(0) {
	TRACE("current node has been initialized to " << current);

}

PseudoDynamic::~PseudoDynamic() {
	// TODO Auto-generated destructor stub
}

void PseudoDynamic::initializeGraph() {
	// do nothing - start with an empty graph
}

void PseudoDynamic::generate() {
	TRACE("calling generate, current node: " << current);
	// if there are no more nodes left from the input graph, stop
	if (this->G->numberOfNodes() == Gstatic.numberOfNodes()) {
		throw std::logic_error("all nodes from the static graph have been generated");
	}


	while (! Gstatic.hasNode(current)) {
		TRACE("static does not have node " << current);
		if (current > Gstatic.upperNodeIdBound()) {
			TRACE(current << " has hit upper id bound for static graph");
			break;
		}
		current += 1;
	}
	node x = this->Gproxy->addNode();
	TRACE("mapping static node " << current << " to dynamic node " << x);
	idmap[current] = x; // map static node u to dynamic node x

	// add edges to already existing neighbors
	Gstatic.forNeighborsOf(current, [&](node v){
		if (v <= current) { // assumption: all nodes below the current one exist in the dynamic graph
			node y = idmap[current];
			node z = idmap[v];
			assert (y != none);
			assert (z != none);
			TRACE("adding edge (" << y << "," << z << ")");
			this->Gproxy->addEdge(y, z);
		}
	});

	// go to next node
	current += 1;
	// trigger time step
	this->Gproxy->timeStep();
}

} /* namespace NetworKit */
