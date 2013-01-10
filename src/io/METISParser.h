/*
 * METISParser.h
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#ifndef METISPARSER_H_
#define METISPARSER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>

#include "../aux/Log.h"


namespace EnsembleClustering {

typedef int64_t node;

class METISParser {

public:

	METISParser();

	virtual ~METISParser();

	/**
	 * Open a METIS graph file.
	 */
	virtual void open(std::string graphPath);

	/**
	 * Get the METIS graph file header
	 */
	virtual std::pair<int, int> getHeader();

	/**
	 * Test if graph file has a next line.
	 */
	virtual bool hasNext();

	/**
	 * Get adjacencies from the next line in the METIS graph file.
	 */
	virtual std::vector<node> getNext();

	/**
	 * Close input file and clean up.
	 */
	virtual void close();

protected:

	std::string graphPath;
	std::ifstream graphFile;
	std::string line;
	int nodeCount;

};

} /* namespace EnsembleClustering */
#endif /* METISPARSER_H_ */
