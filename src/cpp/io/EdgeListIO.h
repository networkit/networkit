/*
 * EdgeListIO.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef EDGELISTIO_H_
#define EDGELISTIO_H_

#include <fstream>
#include <iostream>
#include <string>

#include "GraphReader.h"

namespace NetworKit {

/**
 * @ingroup io
 * A reader for the edge list format used by the LFR benchmark generators, defined as:
 * 		list of edges (nodes are labeled from 1 to the number of nodes;
 * 		the edges are ordered and repeated twice, i.e. source-target and target-source).
 *
 * 	The starting index is a parameter to enable other edge list formats.
 */
class EdgeListIO: public NetworKit::GraphReader {

public:

	/** Default constructor */
	EdgeListIO() = default; //nullary constructor for Python shell

	/**
	 * @param[in]	separator	Character used to separate nodes in an edge line.
	 * @param[in]	firstNode	Index of the first node in the file.
	 */
	EdgeListIO(char separator, node firstNode);

	/**
	 * Given the @a path of an input file, read the graph contained.
	 *
	 * @param[in]	path	Input file path.
	 * @return The graph contained in the file at @a path.
	 */
	Graph read(const std::string& path) override;

		/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _read(std::string& path) {
		return new Graph{std::move(read(path))};
	};

	/**
	 * Write the graph to a file.
	 * @param[in]	G		The graph.
	 * @param[in]	path	The output file path.
	 */
	virtual void write(const Graph& G, std::string path);

protected:

	char separator; 	//!< character separating nodes in an edge line
	node firstNode;
};

} /* namespace NetworKit */
#endif /* EDGELISTIO_H_ */
