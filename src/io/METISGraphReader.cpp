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

	std::cout << "[BEGIN] reading graph G(n=" << n << ", m=" << m << ") from METIS file: " << std::flush;	// status bar follows


	int lc = 0;
	double p = 0.0; // percentage for status bar
	node u = 0;
	while (parser.hasNext()) {
		TRACE("line: " << lc++);
		u += 1;
		std::vector<node> adjacencies = parser.getNext();
		for (node v : adjacencies) {
			if (! G.hasEdge(u, v)) { // edges in METIS file are directed, G edges are undirected
				G.insertEdge(u, v);
			}
		}
		if ((u % 1000) == 0) {
			p = ((double) u / (double) n) * 100;
			std::cout << p << "% ";
		}
	}

	// end status bar
	std::cout << "[DONE]" << std::endl;

	assert (G.numberOfNodes() == n);
	assert (G.numberOfEdges() == m);

	return G;
}

} /* namespace EnsembleClustering */
