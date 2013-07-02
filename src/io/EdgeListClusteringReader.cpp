/*
 * EdgeListClusteringReader.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: forigem
 */

#include "EdgeListClusteringReader.h"

namespace NetworKit {

EdgeListClusteringReader::EdgeListClusteringReader(node firstNode) : firstNode(firstNode) {
	// TODO Auto-generated constructor stub

}

EdgeListClusteringReader::~EdgeListClusteringReader(){
	// TODO Auto-generated destructor stub
}

Clustering EdgeListClusteringReader::read(std::string path) {
	std::ifstream file(path);

	// check if file readable
	if (!file) {
		throw std::runtime_error("invalid clustering file");
	}


	std::vector<cluster> temp;

	// push all cluster ids into vector in order of appearance
	std::string line;
	while(std::getline(file, line)) {
		std::vector<std::string> split = Aux::StringTools::split(line, '\t');
		if (split.size() == 2 && split[0] != "#") {
			cluster c = std::atoi(split[1].c_str());
			temp.push_back(c);
		}
	}

	count n = temp.size();
	Clustering zeta(n);

	#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		zeta[u] = temp[u - firstNode + 1];
	}

	return zeta;
}

} /* namespace NetworKit */
