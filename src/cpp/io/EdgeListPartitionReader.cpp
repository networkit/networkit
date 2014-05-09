/*
 * EdgeListPartitionReader.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: forigem
 */

#include "EdgeListPartitionReader.h"

namespace NetworKit {

EdgeListPartitionReader::EdgeListPartitionReader(node firstNode) : firstNode(firstNode) {
	// TODO Auto-generated constructor stub

}

EdgeListPartitionReader::~EdgeListPartitionReader(){
	// TODO Auto-generated destructor stub
}

Partition EdgeListPartitionReader::read(std::string path) {
	std::ifstream file(path);

	// check if file readable
	if (!file) {
		throw std::runtime_error("invalid clustering file");
	}


	std::vector<index> temp;

	// push all cluster ids into vector in order of appearance
	std::string line;
	while(std::getline(file, line)) {
		std::vector<std::string> split = Aux::StringTools::split(line, '\t');
		if (split.size() == 2 && split[0] != "#") {
			index c = std::atoi(split[1].c_str());
			temp.push_back(c);
		}
	}

	count n = temp.size();
	Partition zeta(n);
	index newOmega = 0;
	#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		newOmega = std::max(newOmega, temp[u - firstNode + 1]);
		zeta[u] = temp[u - firstNode + 1];
	}
	zeta.setUpperBound(newOmega+1);
	return zeta;
}

} /* namespace NetworKit */
