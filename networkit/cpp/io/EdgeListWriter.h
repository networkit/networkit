/*
 * EdgeListWriter.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef EDGELISTWRITER_H_
#define EDGELISTWRITER_H_

#include <fstream>
#include <iostream>
#include <string>

#include "GraphReader.h"

namespace NetworKit {

/**
 * A writer for the edge list format used by the LFR benchmark generators, defined as:
 * 		list of edges (nodes are labelled from 1 to the number of nodes;
 * 		the edges are ordered and repeated twice, i.e. source-target and target-source).
 *
 * 	The starting index is a parameter to enable other edge list formats.
 */
class EdgeListWriter {

public:

	EdgeListWriter() = default; //nullary constructor for Python shell

	/**
	 * @param[in]	separator	character used to separate nodes in an edge line
	 * @param[in]	firstNode	index of the first node in the file
	 */
	EdgeListWriter(char separator, node firstNode);

	/**
	 * Write the graph to a file.
	 * @param[in]	G		the graph
	 * @param[in]	path	the output file path
	 */
	void write(const Graph& G, std::string path);

protected:

	char separator; 	//!< character separating nodes in an edge line
	node firstNode;
};

} /* namespace NetworKit */
#endif /* EDGELISTIO_H_ */
