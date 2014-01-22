/*
 * SNAPEdgeListClusteringReader.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#include "SNAPEdgeListClusteringReader.h"

namespace NetworKit {

SNAPEdgeListClusteringReader::SNAPEdgeListClusteringReader() {
	// TODO Auto-generated constructor stub

}

SNAPEdgeListClusteringReader::~SNAPEdgeListClusteringReader() {
	// TODO Auto-generated destructor stub
}

std::vector<std::set<node>> SNAPEdgeListClusteringReader::read(std::string path) {
	std::ifstream file(path);

	// check if file readable
	if (!file) {
		throw std::runtime_error("invalid clustering file");
	}

	std::vector<std::set<node>> clusters;

	std::string line;
	while(std::getline(file, line)) {
		std::vector<std::string> clusterString = Aux::StringTools::split(line, '\t');
		std::set<node> cluster = {};
		for (std::string nodeString : clusterString) {
			node nodeInCluster = std::atoi(nodeString.c_str());
			cluster.insert(nodeInCluster);
		}
		clusters.push_back(cluster);
	}

	return clusters;
}

} /* namespace NetworKit */
