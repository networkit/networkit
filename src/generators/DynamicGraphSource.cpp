/*
 * DynamicGenerator.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include "DynamicGraphSource.h"

namespace NetworKit {

DynamicGraphSource::DynamicGraphSource() : G(NULL), Gproxy(NULL), graphSet(false), graphInitialized(false) {
	// Graph and GraphEventProxy are set by calling newGraph
}




DynamicGraphSource::~DynamicGraphSource() {
	// TODO Auto-generated destructor stub
}

void DynamicGraphSource::generateWhile(std::function<bool(void)> cont) {
	while (cont()) {
		this->generate();
	}
}

void DynamicGraphSource::generateNodes(count n) {
	auto cont = [&](){
		return (this->G->numberOfNodes() < n);
	};
	this->generateWhile(cont);
}



void DynamicGraphSource::generateEdges(count m) {
	auto cont = [&](){
		return (this->G->numberOfEdges() < m);
	};
	this->generateWhile(cont);
}

GraphEventProxy* DynamicGraphSource::newGraph() {
	this->G = new Graph(0);
	this->Gproxy = new GraphEventProxy(*(this->G));
	// not returning proxy because only generator needs write access to graph
	this->graphSet = true;
	return this->Gproxy;
}

void DynamicGraphSource::generateTimeSteps(count t) {
	auto cont = [&](){
		return (this->G->time() < t);
	};
	this->generateWhile(cont);
}

} /* namespace NetworKit */
