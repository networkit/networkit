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

static inline node fast_string_to_node(std::string::iterator it, const std::string::iterator& end) {
	node val = *it - '0';	// works on ASCII code points
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
	if (parts.size > 2) {
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

	if (flag == 0) {
		while (std::getline(stream, line)) {
			++u;
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
				node v = fast_string_to_node(it1, it2);
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
		}
	} else {
		ERROR("parsing graphs with weight flag " << flag << " not yet supported");
		exit(1);
	}


	return G;
}

} /* namespace NetworKit */
