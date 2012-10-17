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

namespace EnsembleClustering {

METISGraphParser::METISGraphParser() {
	// TODO Auto-generated constructor stub

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
			indices = parseLine(line);
			if (linecount == 0) {
				// handle header line
				// graph data structure
				int n;  // number of nodes
				int m;	// number of edges

				n = indices[0];
				m = indices[1];
				this->initGraph(n, m);
			} else {
				// handle node line
				this->connectNode(currentNode, indices);
				++currentNode;
			}
			++linecount;
		}
	} else {
		std::cout << "unable to open file: " << path << std::endl;
	}

	graphFile.close();

}

void METISGraphParser::initGraph(int n, int m) {

	// TODO: implement

}

void METISGraphParser::connectNode(id v, std::vector<id> indices) {
	std::cout << "Connecting node " << v << " with " << indices.size() << " other nodes" << std::endl;
	// TODO: implement

	int deg = indices.size();


}

} /* namespace EnsembleClustering */
