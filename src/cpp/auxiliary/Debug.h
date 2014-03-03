/*
 * Debug.h
 *
 *  Created on: 23.04.2013
 *      Author: cls
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include <vector>
#include <sstream>
#include <unordered_set>

namespace Aux {


/**
 * String representation of std::vector<T>
 * 		[x, y, ..., z]
 */
template <typename T> std::string vectorToString(const std::vector<T>& vec) {
	std::ostringstream ss;
	ss << '[';
	bool first = true;
	for (const auto& element : vec) {
	    if (!first) {
	        ss << ", ";
	    }
	    ss << element;
	    first = false;
	}
	ss << ']';
	return ss.str();
}

/**
 * String representation of std::unordered_set<T>
 * 		[x, y, ..., z]
 */
template <typename T> std::string setToString(const std::unordered_set<T>& set) {
	std::ostringstream ss;
	ss << '[';
	bool first = true;
	for (const auto& element : set) {
	    if (!first) {
	        ss << ", ";
	    }
	    ss << element;
	    first = false;
	}
	ss << ']';
	return ss.str();
}


} /* namespace Aux */


#endif /* DEBUG_H_ */
