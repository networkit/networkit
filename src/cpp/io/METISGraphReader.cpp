/*
 * METISGraphReader.cpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISGraphReader.h"
#include "METISParser.h"

#include "../auxiliary/Log.h"
#include "../auxiliary/StringTools.h"

namespace NetworKit {

Graph METISGraphReader::read(const std::string& path) {

	METISParser parser(path);

	std::tuple<count, count, index, count> header = parser.getHeader();
	count n = std::get<0>(header);
	count m = std::get<1>(header);
	index fmt = std::get<2>(header);
	count ncon = std::get<3>(header);

	bool weighted;
	if (fmt % 10 == 1) {
		weighted = true;
		DEBUG("graph has been identified as weighted");
	} else {
		weighted = false;
	}
	count ignoreFirst = 0;
	if (fmt / 10 == 1) {
		DEBUG("first ",ncon," value(s) will be ignored");
		ignoreFirst = ncon;
	}

	Graph G(n, weighted);
	std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

	G.setName(graphName);

	INFO("\n[BEGIN] reading graph G(n=", n, ", m=", m, ") from METIS file: ", graphName);	// progress bar follows

	double p = 0.0; // percentage for progress bar
	node u = 0; // begin with 0
	count edgeCounter = 0;
	if (weighted == 0) {
		while (parser.hasNext() && u < n) {
			std::vector<node> adjacencies = parser.getNext(ignoreFirst);
			edgeCounter += adjacencies.size();
			for (index i=0; i < adjacencies.size(); i++) {
				if (adjacencies[i] == 0) {
					ERROR("METIS Node ID should not be 0, edge ignored.");
					continue;
				}
				node v = adjacencies[i] - 1; 	// METIS-indices are 1-based
				assert (v >= 0);
				assert (v < n);
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
	} else {
		while (parser.hasNext() && u < n) {

			std::vector<std::pair<node,double>> adjacencies = parser.getNextWithWeights(ignoreFirst);
			edgeCounter += adjacencies.size();
			DEBUG("node ",u," has ",adjacencies.size()," edges");
			for (index i=0; i < adjacencies.size(); i++) {
				node v = adjacencies[i].first- 1; 	// METIS-indices are 1-based
				double weight = adjacencies[i].second;
				assert (v >= 0);
				if (u <= v && weight > 0) { // self-loops are allowed
					G.addEdge(u, v);
					G.setWeight(u, v, adjacencies[i].second);
					TRACE("(",u,",",v,",",adjacencies[i].second,")");
					assert(adjacencies[i].second > 0);
				}
			}
			u += 1; // next node
			if ((u % ((n + 10)/10)) == 0) {
				p = ((double) (u-1) / (double) n) * 100;
				DEBUG(p, "% ");
			}
		}
	}
	if (G.numberOfEdges() != m) {
		ERROR("METIS file is corrupted: actual number of edges don't match the given number of edges");
	}
	if (edgeCounter / m != 2) {
		WARN("METIS file is corrupted: not every edge is listed twice");
	}

	INFO("\n[DONE]\n");
	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
