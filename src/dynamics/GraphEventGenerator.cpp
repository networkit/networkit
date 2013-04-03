/*
 * GraphEventGenerator.cpp
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#include "GraphEventGenerator.h"

namespace NetworKit {


GraphEventGenerator::GraphEventGenerator(Graph& G) {
	// TODO Auto-generated constructor stub

}

GraphEventGenerator::~GraphEventGenerator() {
	// TODO Auto-generated destructor stub
}

void GraphEventGenerator::generateStream(std::function<bool(void)> done) {

	while (!done()) {


	}
	return;
}

} /* namespace NetworKit */
