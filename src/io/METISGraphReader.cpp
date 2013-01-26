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

	std::cout << "[BEGIN] reading graph G(n=" << n << ", m=" << m << ") from METIS file: " << std::flush;	// progress bar follows


	int lc = 0;
	double p = 0.0; // percentage for progress bar
	node u = 0;
	while (parser.hasNext()) {
		TRACE("line: " << lc++);
		u += 1;
		std::vector<node> adjacencies = parser.getNext();
		for (node v : adjacencies) {
			G.insertEdge(u, v);
		}
		if ((u % 100) == 0) {
			p = ((double) u / (double) n) * 100;
			std::cout << p << "% " << std::flush;
		}
	}

	// end progress bar
	std::cout << "[DONE]" << std::endl;

	assert (G.numberOfNodes() == n);
	assert (G.numberOfEdges() == m);

	return G;
}

} /* namespace EnsembleClustering */
