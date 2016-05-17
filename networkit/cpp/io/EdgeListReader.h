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
#include <map>


#include "GraphReader.h"

namespace NetworKit {

/**
 * @ingroup io
 * A reader for various edge list formats, in which each line contains an edge as
 * two node ids.
 *
 */
class EdgeListReader: public NetworKit::GraphReader {

public:

	EdgeListReader() = default; //nullary constructor for Python shell

	/**
	 * @param[in]	separator	character used to separate nodes in an edge line
	 * @param[in]	firstNode	index of the first node in the file
	 * @param[in]	commentChar	character used to mark comment lines
	 * @param[in]	continuous	boolean to specify, if node ids are continuous
	 * @param[in]	directed	treat graph as directed
	 */
	EdgeListReader(const char separator, const node firstNode, const std::string commentPrefix = "#", const bool continuous = true, const bool directed = false);

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	Graph read(const std::string& path);

	/**
	 * Write the graph to a file.
	 * @param[in]	G		the graph
	 * @param[in]	path	the output file path
	 */
//	virtual void write(const Graph& G, std::string path);

	/**
	 * Return the node map, in case node ids are not continuous
	 */
	std::map<std::string,node> getNodeMap();

protected:
	char separator; 	//!< character separating nodes in an edge line
	std::string commentPrefix;
	node firstNode;
	bool continuous;
	std::map<std::string,node> mapNodeIds;
	bool directed;

private:
	Graph readContinuous(const std::string& path);

	Graph readNonContinuous(const std::string& path);

};

} /* namespace NetworKit */
#endif /* EDGELISTREADER_H_ */
