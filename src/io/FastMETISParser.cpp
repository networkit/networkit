/*
 * FastMETISParser.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISParser.h"

namespace NetworKit {

FastMETISParser::FastMETISParser() {
	// TODO Auto-generated constructor stub

}

FastMETISParser::~FastMETISParser() {
	// TODO Auto-generated destructor stub
}

static inline node fast_string_to_node(std::string::iterator it, const std::string::iterator& end) {
	node val = *it - '0';
	++it;
	while(it != end) {
		val *= 10;
		val += *it  - '0';
		++it;
	}
	return val;
}


std::vector<std::vector<node> > FastMETISParser::parse(std::istream& stream) {
	std::vector<std::vector<node>> data;
	std::string line;
	std::getline(stream, line); // get and discard header
	while(std::getline(stream, line)) {
		std::vector<node> tmp_vec;
		if(line.empty()){
			data.emplace_back(std::move(tmp_vec));
			continue;
		}
		auto it1 = line.begin();
		auto end = line.end();
		auto it2 = std::find(it1, end, ' ');
		if(line.back() == ' '){
			--end;
		}
		while(true) {
			tmp_vec.push_back(fast_string_to_node(it1, it2));
			if(it2 == end) {
				break;
			}
			++it2;
			it1 = it2;
			it2 = std::find(it1, end, ' ');
		}
		data.emplace_back(std::move(tmp_vec));
	}
	return data;
}

} /* namespace NetworKit */
