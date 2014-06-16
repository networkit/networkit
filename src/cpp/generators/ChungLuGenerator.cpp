/*
 * ChungLu.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include "ChungLuGenerator.h"

namespace NetworKit {

ChungLuGenerator::ChungLuGenerator(const std::vector<unsigned int>& degreeSequence) :
		StaticDegreeSequenceGenerator(degreeSequence) {
	sum_deg = std::accumulate(seq.begin(), seq.end(), 0);
	n = (count) seq.size();
}

ChungLuGenerator::~ChungLuGenerator() {

}

Graph ChungLuGenerator::generate() {
	/* Random number in [0, 1] */
	double randVal = 0.0;

	Graph G(n);
	for (index u = 0; u < n; ++u) {
		for (index v = u + 1; v < n; ++v) {
			randVal = Aux::Random::probability();
			/* Probability of edge (u, v): d(u)*d(v)/sum_deg */
			if (randVal < double(seq[u] * seq[v]) / sum_deg) {
				G.addEdge(u, v);
			}
		}
	}

	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
