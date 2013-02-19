/*
 * ParallelAgglomerativeClusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#include "ParallelAgglomerativeClusterer.h"

namespace EnsembleClustering {

ParallelAgglomerativeClusterer::ParallelAgglomerativeClusterer() {
	// TODO Auto-generated constructor stub

}

ParallelAgglomerativeClusterer::~ParallelAgglomerativeClusterer() {
	// TODO Auto-generated destructor stub
}

Clustering ParallelAgglomerativeClusterer::run(Graph& graph) {
	Graph G = graph;
	int attrId = G.addEdgeAttribute_double(0.0);

	do {
		ModularityScoring<double> modScoring(G);
		modScoring.scoreEdges(attrId);

		G.parallelForEdges([&](node u, node v){
			double deltaMod = modScoring.edgeScore(u, v);
			// TODO: where to store edge scores? copy graph?
		});

		// compute matching
		ParallelMatcher parMatcher;
		Matching M = parMatcher.run(G);

		// TODO: contract graph according to matching (and star-like structures)
		MatchingContracter matchingContracter;
		G = matchingContracter.run(G, M);
	} while (false); 	// TODO: repeat until?

	Clustering zeta(G.numberOfNodes());
	// TODO: get clustering induced by contracted graph

	return zeta;
}

} /* namespace EnsembleClustering */
