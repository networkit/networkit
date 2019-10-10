/*
 * LineFileReader.h
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_LINE_FILE_READER_HPP_
#define NETWORKIT_IO_LINE_FILE_READER_HPP_

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
#endif // NETWORKIT_IO_LINE_FILE_READER_HPP_
