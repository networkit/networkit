/*
 * PartitionWriter.h
 *
 *  Created on: 22.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_IO_PARTITION_WRITER_HPP_
#define NETWORKIT_IO_PARTITION_WRITER_HPP_

#include <fstream>

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Write a clustering to a file.
 */
class PartitionWriter {

public:
    virtual void write(Partition& zeta, const std::string& path) const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_PARTITION_WRITER_HPP_
