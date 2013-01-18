/*
 * METISGraphReader.cpp
 *
 *  Created on: 17.01.2013
 *      Author: cls
 */

#include "METISGraphReader.h"

namespace EnsembleClustering {

METISGraphReader::METISGraphReader() {
	// TODO Auto-generated constructor stub

}

METISGraphReader::~METISGraphReader() {
	// TODO Auto-generated destructor stub
}

Graph METISGraphReader::read(std::string path) {

	METISParser parser(path);

	std::pair<int64_t, int64_t> header = parser.getHeader();
	int64_t n = header.first;
	int64_t m = header.second;

	Graph G(n);

	node u = 0;
	while (parser.hasNext()) {
		DEBUG("next line");
		u += 1;
		std::vector<node> adjacencies = parser.getNext();
		for (node v : adjacencies) {
			G.insertEdge(u, v);
		}
	}


	assert (G.numberOfNodes() == n);
	assert (G.numberOfEdges() == m);

	return G;
}

} /* namespace EnsembleClustering */
