/*
 * LineFileReader.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef LINEFILEREADER_H_
#define LINEFILEREADER_H_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace NetworKit {

/**
 * @ingroup io
 * Reads a file and puts each line as a string into a vector
 */
class LineFileReader {
public:
	std::vector<std::string> read(std::string path);
};

} /* namespace NetworKit */
#endif /* LINEFILEREADER_H_ */
