/*
 * METISParser.cpp
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#include "METISParser.h"


namespace EnsembleClustering {


/**
 * Extract a vector of indices from a line in the file.
 *
 * @param[in]	line		line from input file containing node indices
 *
 * @param[out]	indices		node indices extracted from line
 */
static std::vector<node> parseLine(std::string line) {

	std::stringstream stream(line);
	std::string token;
	char delim = ' ';
	std::vector<node> adjacencies;
	node v;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		v = atoi(token.c_str());
		adjacencies.push_back(v);
	}

	return adjacencies;
}


METISParser::METISParser() {

}

METISParser::~METISParser() {
	// TODO Auto-generated destructor stub
}


void METISParser::open(std::string graphPath) {
	TRACE("opening file " << graphPath);
	this->graphPath = graphPath;
	// open METIS graph file
	this->graphFile.open(graphPath.c_str());

	if (graphFile.is_open()) {
		this->nodeCount = 0;
	}  else {
		std::cout << "unable to open file: " << graphPath << std::endl;
	}
}


std::pair<int, int> METISParser::getHeader() {

	// handle header line
	int n;  // number of nodes
	int m;	// number of edges

	std::getline(this->graphFile, this->line);
	std::vector<node> tokens = parseLine(this->line);
	n = tokens[0];
	m = tokens[1];

	TRACE("n = " << n << " m = " << m );

	return std::make_pair(n, m);
}


bool METISParser::hasNext() {
	// if graph file has lines left, return true
	return this->graphFile.good();
}




std::vector<node> METISParser::getNext() {

	bool comment = false;
	do {
		std::getline(this->graphFile, this->line);
		TRACE("reading line: " << line);
		// check for comment line starting with '%'
		if (line[0] == '%') {
			comment = true;
		} else {
			return parseLine(this->line);
		}
	} while (comment);
}


void METISParser::close() {
	// close graph file and clean up
	this->graphFile.close();
	this->line.clear();
}


} /* namespace EnsembleClustering */
