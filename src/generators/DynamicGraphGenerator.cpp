/*
 * DynamicGenerator.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include "DynamicGraphGenerator.h"

namespace NetworKit {

DynamicGraphGenerator::DynamicGraphGenerator(GraphEventProxy& proxy) {
	this->Gproxy = &proxy;
	this->G = this->Gproxy->G;
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

} /* namespace NetworKit */
