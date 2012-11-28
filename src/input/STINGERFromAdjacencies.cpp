/*
 * STINGERFromAdjacencies.cpp
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#include "STINGERFromAdjacencies.h"

#include "log4cxx/logger.h"

extern "C" {
	#include "stinger.h"
}

namespace EnsembleClustering {

STINGERFromAdjacencies::STINGERFromAdjacencies() {
	this->currentNode = 0; // set current node to -1, increment in addAdjacencies
}

STINGERFromAdjacencies::~STINGERFromAdjacencies() {
	// TODO Auto-generated destructor stub
}

void STINGERFromAdjacencies::createGraph() {
	this->G = new Graph(stinger_new());

}

void STINGERFromAdjacencies::addAdjacencies(std::vector<node> adj) {


	this->currentNode++;
	node from = this->currentNode;
	node to;
	for (auto toPtr = adj.begin(); toPtr != adj.end(); ++toPtr) {
		this->G->insertEdge(from, *toPtr);
		LOG4CXX_TRACE(log4cxx::Logger::getRootLogger(), "inserted edge (" << from << "," << *toPtr << ")");
	}

}

stinger* STINGERFromAdjacencies::getSTINGER() {
	return this->G->asSTINGER();
}

Graph* STINGERFromAdjacencies::getGraph() {
	return this->G;
}




} /* namespace EnsembleClustering */
