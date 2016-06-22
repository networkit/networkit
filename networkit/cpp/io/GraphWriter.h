/*
 * GraphWriter.h
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHWRITER_H_
#define GRAPHWRITER_H_

#include "../graph/Graph.h"

#include <fstream>

namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for graph writers.
 */
class GraphWriter {
public:
	virtual ~GraphWriter() = default;

	virtual void write(const Graph& G, const std::string& path) = 0;
};

} /* namespace NetworKit */
#endif /* GRAPHWRITER_H_ */
