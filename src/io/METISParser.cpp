/*
 * METISParser.cpp
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISParser.h"


namespace NetworKit {


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

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		if (token.size() != 0) {
			node v = atoi(token.c_str());
			adjacencies.push_back(v);
		}
	}

	return adjacencies;
}

static std::vector<std::pair<node,double>> parseWeightedLine(std::string line) {

	std::stringstream stream(line);
	std::string token;
	char delim = ' ';
	std::vector<std::pair<node,double>> adjacencies;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		if (token.size() != 0) {
			node v = atoi(token.c_str());
			std::getline(stream, token, delim);
			double weight = atof(token.c_str());
			adjacencies.push_back(std::make_pair(v,weight));
		}

	}

	return adjacencies;
}

static inline std::vector<node> parseLineDIY(std::string line) {
	// TODO: is this faster?

	std::string token;
	char delim = ' ';
	std::vector<node> adjacencies;

	for (char& c : line) {
		if (c == delim) {
			node v = atoi(token.c_str());
			adjacencies.push_back(v);
			token.clear();
		} else if (c == '\n') {
			break;
		} else {
			token.push_back(c);
		}
	}

	return adjacencies;
}


METISParser::METISParser(std::string path) : graphFile(path) {
	if (!(this->graphFile)) {
		ERROR("invalid graph file: " << path);
		throw std::runtime_error("invalid graph file");
	}
}



METISParser::~METISParser() {
	// TODO Auto-generated destructor stub
}


std::tuple<int, int, int> METISParser::getHeader() {

	// handle header line
	int n;  // number of nodes
	int m;	// number of edges
	int weighted; // weighted or unweighted graph

	std::string line = "";
	assert (this->graphFile);

	if (std::getline(this->graphFile, line)) {
		while (line[0] == '%') {
			std::getline(this->graphFile, line);
		}

		std::vector<node> tokens = parseLine(line);
		n = tokens[0];
		m = tokens[1];
		if (tokens.size() == 2) {
			return std::tuple<int, int, int>(n,m,0);
		}
		if (tokens.size() == 3) {
			if (tokens[2] < 2) {
				weighted = tokens[2];
			} else {
				throw std::runtime_error("nodes are weighted");
				return std::tuple<int, int, int>(0,0,0);
			}
		}
		return std::tuple<int, int, int>(n,m,weighted);
	} else {
		ERROR("getline not successful");
		throw std::runtime_error("getting METIS file header failed");
		return std::tuple<int, int, int>(0,0,0);
	}
}




bool METISParser::hasNext() {
	// if graph file has lines left, return true
	return this->graphFile.good();
}


std::vector<node> METISParser::getNext() {

	std::string line;
	bool comment = false;
	do {
		comment = false;
		std::getline(this->graphFile, line);
		// check for comment line starting with '%'
		if (line[0] == '%') {
			comment = true;
		} else {
			return parseLine(line);
		}

	} while (comment);

	throw std::runtime_error("bad METIS file structure");
	std::vector<node> fail;
	return fail;
}

std::vector<std::pair<node,double>> METISParser::getNextWithWeights() {

	std::string line;
	bool comment = false;
	do {
		comment = false;
		std::getline(this->graphFile, line);
		// check for comment line starting with '%'
		if (line[0] == '%') {
			comment = true;
		} else {
			return parseWeightedLine(line);
		}

	} while (comment);

	throw std::runtime_error("bad METIS file structure");
	std::vector<std::pair<node,double>> fail;
	return fail;

}


} /* namespace NetworKit */
