#include "CoverReader.h"

#include <fstream>

NetworKit::Cover NetworKit::CoverReader::read(std::string path, NetworKit::Graph &G)
{
	std::ifstream file;
	file.open(path);
	if (!file.good()) {
		throw std::runtime_error("unable to read from file");
	}
	Cover communities(G.upperNodeIdBound());
	std::string line;
	count i = 0;
	node current;

	while (file.good()) {
		communities.setUpperBound(i+1);
		std::getline(file, line);
		std::stringstream linestream(line);
		while (linestream >> current) {
			communities.addToSubset(i, current);
		}
		++i;
	}

	file.close();
	return communities;
}
