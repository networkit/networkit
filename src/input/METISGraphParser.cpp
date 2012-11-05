/*
 * METISGraphParser.cpp
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#include "METISGraphParser.h"

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "log4cxx/logger.h"



namespace EnsembleClustering {

METISGraphParser::METISGraphParser() {
}

METISGraphParser::~METISGraphParser() {
	// TODO Auto-generated destructor stub
}



/**
 * Extract a vector of indices from a line in the file.
 *
 * @param[in]	line		line from input file containing node indices
 *
 * @param[out]	indices		node indices extracted from line
 */
static std::vector<id> parseLine(std::string line) {


	std::stringstream stream(line);
	std::string token;
	char delim = ' ';
	std::vector<id> indices;
	int id;

	while (std::getline(stream, token, delim)) {
		id = atoi(token.c_str());
		indices.push_back(id);
	}

	return indices;
}


Graph* METISGraphParser::parse(std::string path) {

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "Parsing graph file : " << path);

	std::ifstream graphFile;
	std::string line;
	unsigned int linecount;
	id currentNode;



	graphFile.open(path.c_str());

	if (graphFile.is_open()) {


		linecount = 0;
		currentNode = 0;


		std::vector<id> indices;

		while (graphFile.good()) {
			std::getline(graphFile, line);
			LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "reading line: " << line);
			// check for comment line starting with '%'
			if (line[0] == '%') {
				// omit comment line
			} else {
				LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "linecount: " << linecount);
				indices = parseLine(line);
				if (linecount == 0) {
					// handle header line
					// graph data structure
					int n;  // number of nodes
					int m;	// number of edges

					n = indices[0];
					m = indices[1];

					LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "n = " << n << " m = " << m )

					this->initGraph(n, m);

				} else {
					// handle node line
					// FIXME:adapt to new graph data structure
					// this->graphData->connectNode(currentNode, indices);
					++currentNode;
				}
				++linecount;
			}

		} // end while

		LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "finished reading file");


	} else {
		std::cout << "unable to open file: " << path << std::endl;
	}

	graphFile.close();


	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "Done: Graph file parsed");

	// TODO: assemble graph
	// Graph* G = new Graph();
	// return G;

}

void METISGraphParser::initGraph(int n, int m) {

	// FIXME:
	// this->graphData = new EdgeTripleGraph(n, m);


	// FIXME:
	// LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "created new graph data structure with n=" << this->graphData->n << " and m=" << this->graphData->m );

}



} /* namespace EnsembleClustering */
