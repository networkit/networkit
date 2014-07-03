/*
 * GraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHREADER_H_
#define GRAPHREADER_H_

#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"


namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for graph readers.
 */
class GraphReader {
public:
	virtual ~GraphReader() = default;

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	virtual Graph read(const std::string& path) = 0;

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _read(std::string& path) {
		return new Graph{std::move(read(path))};
	};
};

} /* namespace NetworKit */
#endif /* GRAPHREADER_H_ */
