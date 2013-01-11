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

Clustering ParallelAgglomerativeClusterer::run(Graph& G) {

	do {
		// TODO: score edges with deltaMod
		ModularityScoring modScoring(G);	// TODO: make this EdgeScoring to allow other scoring measures

		G.forallEdges([&](node u, node v){
			double deltaMod = modScoring.scoreEdge(u, v);
			// TODO: where to store edge scores? copy graph?
		}, "parallel", "readonly");

		// compute matching
		ParallelMatcher parMatcher;
		Matching M = parMatcher.run(G);

		// TODO: contract graph according to matching (and star-like structures)
		MatchingContracter matchingContracter;
		Graph Gcon = matchingContracter.run(G, M);

	} while (false); 	// TODO: repeat until?

	Clustering zeta(G.numberOfNodes());
	// TODO: get clustering induced by contracted graph

	return zeta;
}

} /* namespace EnsembleClustering */
