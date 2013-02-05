/*
 * GraphFromAdjacencies.cpp
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#include "GraphFromAdjacencies.h"



namespace EnsembleClustering {

GraphFromAdjacencies::GraphFromAdjacencies() {
	this->currentNode = 1; //<! first node index;  METIS graph files have 1-based indices
}

GraphFromAdjacencies::~GraphFromAdjacencies() {
	// TODO Auto-generated destructor stub
}

void GraphFromAdjacencies::readHeader(std::pair<int64_t, int64_t> header) {
	this->G.extendNodeRange(header.first);	// allocate node range to number of nodes
}

void GraphFromAdjacencies::addAdjacencies(std::vector<node> adj) {


	node from = this->currentNode;
	for (auto toPtr = adj.begin(); toPtr != adj.end(); ++toPtr) {
		this->G.insertEdge(from, *toPtr);
		TRACE("inserted edge (" << from << "," << *toPtr << ")");
	}
	this->currentNode++;
}



Graph GraphFromAdjacencies::getGraph() {
	return this->G;
}




} /* namespace EnsembleClustering */
