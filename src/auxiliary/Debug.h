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

namespace Aux {


/**
 * Get string representation of vector.
 */
template <typename T> std::string vectorToString(std::vector<T>& vec) {
	std::stringstream out;
	out << "[";
	for (T element : vec) {
		out << element << ",";
	}
	out << "]";
	return out.str();
}


} /* namespace Aux */


#endif /* DEBUG_H_ */
