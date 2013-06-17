/*
 * DynamicGenerator.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include "DynamicGraphGenerator.h"

namespace NetworKit {

DynamicGraphGenerator::DynamicGraphGenerator() : G(NULL), Gproxy(NULL), graphSet(false), graphInitialized(false) {
	// Graph and GraphEventProxy are set by calling newGraph
}




DynamicGraphGenerator::~DynamicGraphGenerator() {
	// TODO Auto-generated destructor stub
}

void DynamicGraphGenerator::generateWhile(std::function<bool(void)> cont) {
	while (cont()) {
		this->generate();
	}
}

void DynamicGraphGenerator::generateNodes(count n) {
	auto cont = [&](){
		return (this->G->numberOfNodes() < n);
	};
	this->generateWhile(cont);
}



void DynamicGraphGenerator::generateEdges(count m) {
	auto cont = [&](){
		return (this->G->numberOfEdges() < m);
	};
	this->generateWhile(cont);
}

GraphEventProxy* DynamicGraphGenerator::newGraph() {
	this->G = new Graph(0);
	this->Gproxy = new GraphEventProxy(*(this->G));
	// not returning proxy because only generator needs write access to graph
	this->graphSet = true;
	return this->Gproxy;
}

void DynamicGraphGenerator::generateTimeSteps(count t) {
	auto cont = [&](){
		return (this->G->time() < t);
	};
	this->generateWhile(cont);
}

} /* namespace NetworKit */
