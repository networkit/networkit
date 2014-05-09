/*
 * PartitionWriter.cpp
 *
 *  Created on: 22.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "PartitionWriter.h"

namespace NetworKit {

PartitionWriter::PartitionWriter() {
	// TODO Auto-generated constructor stub

}

PartitionWriter::~PartitionWriter() {
	// TODO Auto-generated destructor stub
}

void PartitionWriter::write(Partition& zeta, std::string path) const {
	std::ofstream file;
	file.open(path.c_str());

	zeta.forEntries([&](node v, index c){
		file << c << std::endl;
	});

	file.close();
}

} /* namespace NetworKit */
