/*
 * KONECTReader.h
 *
 */

#ifndef KONECTGRAPHREADER_H_
#define KONECTGRAPHREADER_H_

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>


#include "GraphReader.h"

namespace NetworKit {

/*
 * KONECTGraphReader.cpp
 * 
 * Reader for the KONECT graph format, 
 * based on the EdgeListReader.
 * 
 * The KONECT format is described in detail in 
 * http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
 */
class KONECTGraphReader: public NetworKit::GraphReader {

public:

	KONECTGraphReader() = default; //nullary constructor for Python shell

	/**
	 * @param[in]	ignoreLoops	ignores loops in the input graph file, if set to true
	 * @param[in]	separator	character used to separate values of a line
	 */
	KONECTGraphReader(char separator, bool ignoreLoops=false);

	/**
	 * Given the path of an input file, read the graph contained.
	 *
	 * @param[in]	path	input file path
	 */
	Graph read(const std::string& path);


protected:
	char separator; 	//!< character separating nodes in an edge line
	std::string commentPrefix;
	node firstNode;
	bool continuous;
//	std::unordered_map<index,node> mapNodeIds;
	bool ignoreLoops;

private:
	Graph readContinuous(const std::string& path);

};

} /* namespace NetworKit */
#endif /* KONECTGRAPHREADER_H_ */
