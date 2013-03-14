/*
 * GraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHREADER_H_
#define GRAPHREADER_H_

#include "../graph/Graph.h"

namespace EnsembleClustering {

class GraphReader {

public:

	GraphReader();

	virtual ~GraphReader();

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	virtual Graph read(std::string path) = 0;
};

} /* namespace EnsembleClustering */
#endif /* GRAPHREADER_H_ */
