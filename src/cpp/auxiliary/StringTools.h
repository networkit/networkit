/*
 * StringTools.h
 *
 *  Created on: 11.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef STRINGTOOLS_H_
#define STRINGTOOLS_H_

#include <vector>
#include <string>
#include <sstream>

namespace Aux {

/**
 * Missing string functions.
 */
namespace StringTools {

/**
 * Split a string at delimiter and return vector of parts.
 *
 */
inline std::vector<std::string> split(const std::string& s, char delim = ' ') {
	std::stringstream stream(s);
	std::string token;
	std::vector<std::string> tokens;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		tokens.push_back(token);
	}

	return tokens;
}

} /* namespace StringTools */
} /* namespace Aux */
#endif /* STRINGTOOLS_H_ */
