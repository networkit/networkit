/*
 * LineFileReader.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "LineFileReader.h"

namespace NetworKit {
std::vector<std::string> LineFileReader::read(std::string path) {
	std::ifstream file;
	std::string line; // the current line
	file.open(path);

	std::vector<std::string> data;

	while (file.good()) {
		std::getline(file, line);
		data.push_back(line);
	}

	file.close();

	return data;
}

} /* namespace NetworKit */
