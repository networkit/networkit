/*
 * FastMETISParser.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISParserDouble.h"

#include <fstream>

namespace NetworKit {

FastMETISParserDouble::FastMETISParserDouble() {
	// TODO Auto-generated constructor stub

}

FastMETISParserDouble::~FastMETISParserDouble() {
	// TODO Auto-generated destructor stub
}
#if 1
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


NetworKit::Graph FastMETISParserDouble::parse(const std::string& path) {
	std::ifstream stream(path);
	std::string line;
	std::getline(stream, line); // get header
	count n;
	count m;
	int flag; // weight flag
	std::tie(n, m, flag) = parseHeader(line);
	if (flag > 1) return Graph(0); // return empty graph in case of weighted nodes
	bool weighted = flag % 10;
	Graph G(n, weighted);

	std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

	G.setName(graphName);

	node u = 0;

	// XXX: if the first character of a line is the % char, comment lines will be ignored

	std::string current;
	node v = 0;
	if (!weighted) {
		DEBUG("reading unweighted graph");
		// unweighted edges
		while (std::getline(stream, line)) {
			if (line.empty() || (!line.empty() && line[0] == '%') ) {
				continue;
			}
			// determine index of last character of last edge weight
			count end = line.length()-1;
			while (line[end] == ' ' && end > 0) {
				--end;
			}

			for(count i = 0; i <= end; ++i) {
				if (line[i] != ' ') {
					current += line[i];
				}
				if (!current.empty() && (line[i] == ' ' || i == end)) {
					v = std::stoi(current);
					current = "";
					if (u < v) G.addEdge(u,v-1);
				}
			}
			++u;
		}
	} else {
		DEBUG("reading weighted graph");
		double weight = 0.0;
		bool decider = false;
		// weighted edges - WARNING: supports only non-negative integer weights
		while (std::getline(stream, line)) {
			if (line.empty() || (!line.empty() && line[0] == '%') ) {
				continue;
			}
			// determine index of last character of last edge weight
			count end = line.length()-1;
			//DEBUG("line length: ", end+1, " with char: ", line[end]);
			while (line[end] == ' ' && end > 0) {
				--end;
			}
			//DEBUG("set end to: ", end, " with char: ", line[end]);
			for(count i = 0;  i <= end; ++i) {
				if (line[i] != ' ') {
					//DEBUG("adding ",line[i], " to current");
					current += line[i];
				}
				if (!current.empty() && (line[i] == ' ' || i == end)) {
				//} else if (!current.empty()) {
					if (!decider) {
						//DEBUG("std::stoi(", current, ")");
						v = std::stoi(current);
					} else {
						//DEBUG("std::stod(", current, ")");
						weight = std::stod(current);
						if (u < v) G.addEdge(u,v-1, weight);
						//DEBUG("edge to graph: ", u, ", ", v-1, ", ", weight);
					}
					decider = !decider;
					current = "";
				}
			}
			++u;
		}
	}


	return G;
}

#else 
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

static inline double fast_string_to_double(std::string::iterator it, const std::string::iterator& end) {
	double val = 0.0;
	bool negative = false;
	if (*it == '-') {
		negative = true;
		++it;
	}
	while (*it >= '0' && *it <= '9' && it != end) { //it != end
		val = (val*10.0) + (*it - '0');
		++it;
	}
	if (*it == '.' && it != end) {
		double afterComma = 0.0;
		int exp = 0;
		++it;
		while (*it >= '0' && *it <= '9' && it != end) { //it != end
			afterComma = (afterComma*10.0) + (*it - '0');
			++exp;
			++it;
		}
		val += afterComma / std::pow(10.0, exp);
	}
	if (negative) {
		val = -val;
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


NetworKit::Graph FastMETISParserDouble::parse(const std::string& path) {
	std::ifstream stream(path);
	std::string line;
	std::getline(stream, line); // get header
	count n;
	count m;
	int flag; // weight flag
	std::tie(n, m, flag) = parseHeader(line);
	bool weighted = flag % 10;
	Graph G(n, weighted);
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
			std::stringstream linestream(line);
			node v;
			double weight;
			while (linestream >> v >> weight) {
				if (u < v) {
					G.addEdge(u,v-1, (edgeweight) weight);
				}
			}
/*			auto it1 = line.begin();
			auto end = line.end();
			auto it2 = std::find(it1, end, ' ');
			if (line.back() == ' ') { // if line ends with one white space, do this trick...
				--end;
			}
			while (true) {
				TRACE("try to read new node");
				node v = (fast_string_to_integer(it1, it2) - 1);
				TRACE("read node " , v);
				// advance
				++it2;
				it1 = it2;
				it2 = std::find(it1, end, ' ');
				TRACE("char at it1: ",*it1," and at it2: ",*it2," and their difference: ",&*it2-&*it1);
				// get weight
				double weight = (fast_string_to_double(it1, it2));
				//double weight = std::stod(line.sub);
				TRACE("read weight " , weight);
				++it2;
				it1 = it2;
				TRACE("find new it2");
				it2 = std::find(it1, end, ' ');
				TRACE("found new it2");
				TRACE("char at it1: ",*it1," and at it2: ",*it2," and their difference: ",&*it2-&*it1);
			
				if (u < v) {
					G.addEdge(u, v, (edgeweight) weight);
				}
				TRACE("new node added");
				if (it2 == end) {
					break;
				}
			}*/
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
#endif

} /* namespace NetworKit */
