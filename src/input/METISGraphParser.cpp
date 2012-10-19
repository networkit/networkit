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


#include "../graph/EdgeTripleGraphData.h"

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


Graph METISGraphParser::parse(std::string path) {

	std::ifstream graphFile;
	std::string line;
	unsigned int linecount;
	id currentNode;



	graphFile.open(path.c_str());

	if (graphFile.is_open()) {


		linecount = 0;
		currentNode = 0;
		std::vector<id> indices; // the integer indices for the current line

		while (graphFile.good()) {
			std::getline(graphFile, line);
			// check for comment line starting with '%'
			if (line[0] == '%') {
				// omit comment line
			} else {
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
					this->graphData->connectNode(currentNode, indices);
					++currentNode;
				}
				++linecount;
			}

		}
	} else {
		std::cout << "unable to open file: " << path << std::endl;
	}

	graphFile.close();

	// TODO: return graph

}

void METISGraphParser::initGraph(int n, int m) {

	this->graphData = new EdgeTripleGraphData(n, m);


	LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "created new graph data structure with n=" << this->graphData->n << " and m=" << this->graphData->m );

}



} /* namespace EnsembleClustering */
