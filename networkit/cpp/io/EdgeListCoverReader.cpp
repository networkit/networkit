#include "EdgeListCoverReader.h"

#include <fstream>
#include "../graph/Graph.h"
#include "../auxiliary/StringTools.h"


NetworKit::EdgeListCoverReader::EdgeListCoverReader(NetworKit::node firstNode) : firstNode(firstNode) {}

NetworKit::Cover NetworKit::EdgeListCoverReader::read(std::string path, Graph &G) {
	std::ifstream file;
	file.open(path);
	if (!file.good()) {
		throw std::runtime_error("unable to read from file");
	}
	Cover communities(G.upperNodeIdBound());

	std::string line;
	index omega = 0;

	while(std::getline(file, line)) {
		if (line.substr(0, 1) != "#") {
			std::stringstream linestream(line);
			index c, v;
			linestream >> v;
			// NetworKit uses zero-based node ids, adapt input accordingly
			v -= firstNode;
			if (v >= G.upperNodeIdBound()) {
				throw std::runtime_error("Node id in the input too large");
			}

			while (linestream >> c) {
				if (c > omega) {
					communities.setUpperBound(c+1);
					omega = c;
				}

				communities.addToSubset(c, v);
			}
		}
	}

	file.close();
	return communities;

}
