/*
 * BinaryPartitionReader.h
 *
 *  Created on: 12.04.2017
 *      Author: Michael Hamann <michael.hamann@kit.edu>
 */

#ifndef BINARYPARTITIONREADER_H_
#define BINARYPARTITIONREADER_H_

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class BinaryPartitionReader {

public:
	/**
	 * Construct a binary partition reader.
	 *
	 * @param[in]	width	The integer width. Supported values: 4 and 8.
	 */
	BinaryPartitionReader(uint8_t width = 4);
	

	/**
	 * Read a partition from a file. File format:
	 * 		list of (unsigned) integer partition ids, one for every node
	 *
	 * @param[in]	path	Path to file.
	 */
	virtual Partition read(const std::string& path);
private:
	uint8_t width;
};

} /* namespace NetworKit */
#endif /* BINARYPARTITIONREADER_H_ */
