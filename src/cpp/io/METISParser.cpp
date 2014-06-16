/*
 * METISParser.cpp
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISParser.h"

#include "../auxiliary/Enforce.h"

namespace NetworKit {


/**
 * Extract a vector of indices from a line in the file.
 *
 * @param[in]	line		line from input file containing node indices
 *
 * @param[out]	indices		node indices extracted from line
 */
static std::vector<node> parseLine(const std::string& line) {

	std::stringstream stream(line);
	std::string token;
	char delim = ' ';
	std::vector<uint64_t> adjacencies;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		if (token.size() != 0) {
			node v = stoull(token);
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
			node v = stoull(token);
			std::getline(stream, token, delim);
			double weight = stod(token);
			adjacencies.push_back(std::make_pair(v,weight));
		}

	}

	return adjacencies;
}


METISParser::METISParser(std::string path) : graphFile(path) {
	if (!(this->graphFile)) {
		ERROR("invalid graph file: " , path);
		throw std::runtime_error("invalid graph file");
	}
}



METISParser::~METISParser() {
	// TODO Auto-generated destructor stub
}


std::tuple<count, count, index> METISParser::getHeader() {

	// handle header line
	count n;  // number of nodes
	count m;	// number of edges
	index weighted; // weighted or unweighted graph

	std::string line = "";
	Aux::enforce (this->graphFile);

	if (std::getline(this->graphFile, line)) {
		while (line[0] == '%') {
			std::getline(this->graphFile, line);
		}

		std::vector<uint64_t> tokens = parseLine(line);
		n = tokens[0];
		m = tokens[1];
		if (tokens.size() == 2) {
			return std::tuple<count, count, index>(n,m,0);
		}
		if (tokens.size() == 3) {
			if (tokens[2] < 2) {
				weighted = tokens[2];
			} else {
				throw std::runtime_error("nodes are weighted");
				return std::tuple<count, count, index>(0,0,0);
			}
		}
		return std::tuple<count, count, index>(n,m,weighted);
	} else {
		ERROR("getline not successful");
		throw std::runtime_error("getting METIS file header failed");
		return std::tuple<count, count, index>(0,0,0);
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
