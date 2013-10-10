/*
 * FastMETISGraphReader.cpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#include "FastMETISGraphReader.h"

namespace NetworKit {

FastMETISGraphReader::FastMETISGraphReader() {
	// TODO Auto-generated constructor stub

}

FastMETISGraphReader::~FastMETISGraphReader() {
	// TODO Auto-generated destructor stub
}

Graph FastMETISGraphReader::read(std::string path) {
	FastMETISParser parser;
	std::ifstream file;
	file.open(path);
	std::string header;
	std::getline(file, header);
	std::cout << header << std::endl;
	count n = std::stoi(Aux::StringTools::split(header, ' ')[0]);
	count m = std::stoi(Aux::StringTools::split(header, ' ')[1]);
	file.close();
	file.open(path);
	std::vector<std::vector<node> > adja = parser.parse(file);

	Graph G(n);
	G.forNodes([&](node u) {
		for (node v : adja[u]) {
			if (u < v) {
				G.addEdge(u, v);
			}
		}
	});
	return G;
}

} /* namespace NetworKit */
