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


} /* namespace NetworKit */
