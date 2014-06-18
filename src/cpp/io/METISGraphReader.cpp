/*
 * METISGraphReader.cpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISGraphReader.h"

namespace NetworKit {

Graph METISGraphReader::read(const std::string& path) {

	METISParser parser(path);

	std::tuple<count, count, index> header = parser.getHeader();
	count n = std::get<0>(header);
	count m = std::get<1>(header);
	index weightMarker = std::get<2>(header);
	//std::tie(n,m,weighted) = header;
	// TODO: std::tie(n, m, weighted) = header

	bool weighted;
	if (weightMarker == 0) {
		weighted = false;
	} else {
		weighted = true;
	}
	Graph G(n, weighted);

	std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

	G.setName(graphName);

	INFO("\n[BEGIN] reading graph G(n=", n, ", m=", m, ") from METIS file: ", graphName);	// progress bar follows

	double p = 0.0; // percentage for progress bar
	node u = 0; // begin with 0

	if (weighted == 0) {
		while (parser.hasNext()) {
			std::vector<node> adjacencies = parser.getNext();
			for (index i=0; i < adjacencies.size(); i++) {

				node v = adjacencies[i] - 1; 	// METIS-indices are 1-based
				assert (v >= 0);
				if (u <= v) { // self-loops are allowed
					G.addEdge(u, v);
				}
			}
			u++; // next node
			if ((u % ((n + 10)/10)) == 0) {
				p = ((double) (u-1) / (double) n) * 100;
				DEBUG(p, "% ");
			}
		}
		INFO("\n[DONE]\n");
		return G;
	} else {
		while (parser.hasNext()) {

			std::vector<std::pair<node,double>> adjacencies = parser.getNextWithWeights();

			for (index i=0; i < adjacencies.size(); i++) {
				node v = adjacencies[i].first- 1; 	// METIS-indices are 1-based
				assert (v >= 0);
				if (u <= v) { // self-loops are allowed
					G.addEdge(u, v);
					G.setWeight(u, v, adjacencies[i].second);
					assert(adjacencies[i].second > 0);
				}
			}
			u += 1; // next node
			if ((u % ((n + 10)/10)) == 0) {
				p = ((double) (u-1) / (double) n) * 100;
				DEBUG(p, "% ");
			}
		}
		INFO("\n[DONE]\n");
		return G;
	}
}

} /* namespace NetworKit */
