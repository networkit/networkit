/*
 * PartitionReader.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "PartitionReader.h"

namespace NetworKit {

Partition PartitionReader::read(std::string path) {

	std::ifstream file(path);

	// check if file readable
	if (!file) {
		throw std::runtime_error("invalid clustering file");
	}


	std::vector<index> temp;

	// push all cluster ids into vector in order of appearance
	std::string line;
	while(std::getline(file, line)) {
		index c = std::atoi(line.c_str());
		temp.push_back(c);
	}

	count n = temp.size();
	Partition zeta(n);

	#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		zeta[u] = temp[u];
	}

	return zeta;
}

} /* namespace NetworKit */
