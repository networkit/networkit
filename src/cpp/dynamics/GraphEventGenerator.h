/*
 * GraphEventGenerator.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTGENERATOR_H_
#define GRAPHEVENTGENERATOR_H_

#include <functional>

#include "GraphEvent.h"

namespace NetworKit {

class GraphEventGenerator {

public:

	GraphEventGenerator(Graph& G);

	virtual ~GraphEventGenerator();

	virtual void generateStream(std::function<bool(void)> done);
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTGENERATOR_H_ */
