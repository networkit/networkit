/*
 * PartitionReader.h
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_IO_PARTITION_READER_HPP_
#define NETWORKIT_IO_PARTITION_READER_HPP_

#include <fstream>

#include <networkit/structures/Partition.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class PartitionReader {

public:
    /**
     * Read a clustering from a file. File format:
     * 		line n contains cluster id of node (n - 1)
     *
     * @param[in]	path	Path to file.
     */
    virtual Partition read(std::string path);
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_PARTITION_READER_HPP_
