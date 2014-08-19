/*
 * PartitionWriter.h
 *
 *  Created on: 22.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef PARTITIONWRITER_H_
#define PARTITIONWRITER_H_

#include <fstream>

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Write a clustering to a file.
 */
class PartitionWriter {

public:

	PartitionWriter() = default;

	virtual ~PartitionWriter() = default;

	virtual void write(Partition& zeta, const std::string& path) const;
};

} /* namespace NetworKit */
#endif /* CLUSTERINGWRITER_H_ */
