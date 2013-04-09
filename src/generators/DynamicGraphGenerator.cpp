/*
 * DynamicGenerator.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include "DynamicGraphGenerator.h"

namespace NetworKit {

DynamicGraphGenerator::DynamicGraphGenerator() {
	// TODO Auto-generated constructor stub

}

DynamicGraphGenerator::~DynamicGraphGenerator() {
	// TODO Auto-generated destructor stub
}

void DynamicGraphGenerator::setProxy(GraphEventProxy& proxy) {
	this->proxy = &proxy;
	this->G = this->proxy->G;
}

} /* namespace NetworKit */
