/*
 * PartitionWriter.cpp
 *
 *  Created on: 22.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "PartitionWriter.h"

namespace NetworKit {

void PartitionWriter::write(Partition& zeta, const std::string& path) const {
	std::ofstream file{path};

	zeta.forEntries([&](node, index c){
		file << c << '\n';
	});
}

} /* namespace NetworKit */
