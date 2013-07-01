/*
 * EdgeListReader.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef EDGELISTREADER_H_
#define EDGELISTREADER_H_

#include <fstream>
#include <iostream>
#include <string>

#include "GraphReader.h"

namespace NetworKit {

/**
 * A reader for the edge list format used by the LFR benchmark generators, defined as:
 * 		list of edges (nodes are labelled from 1 to the number of nodes;
 * 		the edges are ordered and repeated twice, i.e. source-target and target-source).
 *
 * 	The starting index is a parameter to enable other edge list formats.
 */
class EdgeListReader: public NetworKit::GraphReader {

public:

	/**
	 * @param[in]	firstNode	index of the first node in the file
	 */
	EdgeListReader(node firstNode = 1);

	virtual ~EdgeListReader();

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	virtual Graph read(std::string path);

protected:

	node firstNode;
};

} /* namespace NetworKit */
#endif /* EDGELISTREADER_H_ */
