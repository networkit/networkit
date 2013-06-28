/*
 * DynamicDGSParser.h
 *
 *  Created on: Jun 17, 2013
 *      Author: forigem
 */

#ifndef DYNAMICDGSPARSER_H_
#define DYNAMICDGSPARSER_H_

#include <fstream>
#include <vector>


#include "DynamicGraphSource.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

class DynamicDGSParser: public NetworKit::DynamicGraphSource {
public:
	DynamicDGSParser(std::string path);
	virtual ~DynamicDGSParser();

	/**
	 * The generator may expect the graph to be in a certain initial state. Call this method first.
	 */
	virtual void initializeGraph();


	/**
	 * Perform one generative step - as defined by the implementation.
	 */
	virtual void generate();

protected:
	bool graphInitialized;	//!< true if initializeGraph has been called and graph has been properly initialized
	std::unordered_map<std::string, node> nodeNames;
	std::vector<std::vector<std::string>> nodeCategories;
	std::vector<std::string> nodeDates;
	std::ifstream dgsFile;

};

} /* namespace NetworKit */
#endif /* DYNAMICDGSPARSER_H_ */
