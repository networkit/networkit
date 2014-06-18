/*
 * GraphEventGenerator.cpp
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#include "GraphEventGenerator.h"

namespace NetworKit {


GraphEventGenerator::GraphEventGenerator(Graph& G) {
}

GraphEventGenerator::~GraphEventGenerator() {
}

void GraphEventGenerator::generateStream(std::function<bool(void)> done) {

	while (!done()) {


	}
	return;
}

} /* namespace NetworKit */
