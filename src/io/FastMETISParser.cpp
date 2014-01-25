/*
 * FastMETISParser.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISParser.h"

#include <fstream>

namespace NetworKit {

FastMETISParser::FastMETISParser() {
	// TODO Auto-generated constructor stub

}

FastMETISParser::~FastMETISParser() {
	// TODO Auto-generated destructor stub
}

static inline uint64_t fast_string_to_integer(std::string::iterator it, const std::string::iterator& end) {
	uint64_t val = *it - '0';	// works on ASCII code points
	++it;
	while (it != end) {
		val *= 10;
		val += *it  - '0';
		++it;
	}
	return val;
}

static inline std::tuple<count, count, int> parseHeader(const std::string& header) {
	count n;
	count m;
	int flag;

	std::vector<std::string> parts = Aux::StringTools::split(header);
	n = std::stoi(parts[0]);
	m = std::stoi(parts[1]);
	if (parts.size() > 2) {
		flag = std::stoi(parts[2]);
	} else {
		flag = 0;
	}


	return std::make_tuple(n, m, flag);
}


NetworKit::Graph FastMETISParser::parse(const std::string& path) {
	std::ifstream stream(path);
	std::string line;
	std::getline(stream, line); // get header
	count n;
	count m;
	int flag; // weight flag
	std::tie(n, m, flag) = parseHeader(line);
	Graph G(n);
	node u = 0;

	// TODO: handle comment lines

	if (flag == 0) {
		DEBUG("reading unweighted graph");
		// unweighted edges
		while (std::getline(stream, line)) {
			if (line.empty()) {
				continue;
			}
			auto it1 = line.begin();
			auto end = line.end();
			auto it2 = std::find(it1, end, ' ');
			if (line.back() == ' ') { // if line ends with one white space, do this trick...
				--end;
			}
			while (true) {
				node v = (fast_string_to_integer(it1, it2) - 1);
				if (u < v) {
					G.addEdge(u, v);
				}
				if (it2 == end) {
					break;
				}
				++it2;
				it1 = it2;
				it2 = std::find(it1, end, ' ');
			}
			++u; // next node
		}
	} else if (flag == 1) {
		DEBUG("reading weighted graph");
		// weighted edges - WARNING: supports only non-negative integer weights
		while (std::getline(stream, line)) {
			if (line.empty()) {
				continue;
			}
			auto it1 = line.begin();
			auto end = line.end();
			auto it2 = std::find(it1, end, ' ');
			if (line.back() == ' ') { // if line ends with one white space, do this trick...
				--end;
			}
			while (true) {
				node v = (fast_string_to_integer(it1, it2) - 1);
				TRACE("read node " , v);
				// advance
				++it2;
				it1 = it2;
				it2 = std::find(it1, end, ' ');
				// get weight
				uint64_t weight = (fast_string_to_integer(it1, it2));
				TRACE("read weight " , weight);
				++it2;
				it1 = it2;
				it2 = std::find(it1, end, ' ');
				if (u < v) {
					G.addEdge(u, v, (edgeweight) weight);
				}
				if (it2 == end) {
					break;
				}
			}
			++u; // next node
		}
		
	} else if (flag == 11) {
		ERROR("weighted nodes are not supported");
		throw std::runtime_error("weighted nodes are not supported");
	} else {
		ERROR("invalid weight flag in header of METIS file: " , flag);
		throw std::runtime_error("invalid weight flag in header of METIS file");
	}


	return G;
}

} /* namespace NetworKit */
