/*
 * METISParser.h
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef METISPARSER_H_
#define METISPARSER_H_

#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include "../Globals.h"

namespace NetworKit {

/**
 * @ingroup io
 * Parser for the METIS file format.
 */
class METISParser {

protected:

	std::ifstream graphFile;


public:
	/**
	 * Constructor for the METIS Parser.
	 *
	 * @param[in]	path	file path 
	 */
	METISParser(std::string path);


	/**
	 * Get the METIS graph file header
	 */
	std::tuple<count, count, index, count> getHeader();

	/**
	 * Test if graph file has a next line.
	 */
	bool hasNext();

	/**
	 * Get adjacencies from the next line in the METIS graph file.
	 *
	 * @param[in]	ignoreFirst	number of values to ignore [in case the METIS file contains node weightes]
	 */
	std::vector<node> getNext(count ignoreFirst = 0);

	/**
	 * Get adjacencies with edge weights from the next line in the METIS graph file.
	 *
	 * @param[in]	ignoreFirst	number of values to ignore [in case the METIS file contains node weightes]
	 */
	std::vector<std::pair<node,double>> getNextWithWeights(count ignoreFirst = 0);

};
} /* namespace NetworKit */
#endif /* METISPARSER_H_ */
