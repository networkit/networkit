/*
 * ParallelAgglomerativeClusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu),
 *      		Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include "ParallelAgglomerativeClusterer.h"
#include "../scoring/ModularityScoring.h"
#include "../matching/PathGrowingMatcher.h"
#include "../coarsening/MatchingContracter.h"
#include "../coarsening/ClusteringProjector.h"

namespace NetworKit {

Partition ParallelAgglomerativeClusterer::run(const Graph& graph) {
	// copy graph because we make changes due to merges
	Graph G(graph.numberOfNodes(), true); // make weighted copy
	graph.forWeightedEdges([&](node u, node v, edgeweight w){
		G.addEdge(u, v, w);
	});

	std::vector<std::vector<node> > mapHierarchy;

	bool repeat = true;
	do {
		// prepare attributes for scoring
		int attrId = G.addEdgeAttribute_double(0.0);

		// perform scoring
		TRACE("before scoring graph of size " , G.numberOfNodes());
		ModularityScoring<double> modScoring(G);
		modScoring.scoreEdges(attrId);

		// FIXME: so far only sequential
		// compute matching
		PathGrowingMatcher parMatcher;
		Matching M = parMatcher.run(G);

		// contract graph according to matching, TODO: (and star-like structures)
		MatchingContracter matchingContracter;
		auto GandMap = matchingContracter.run(G, M);

		// determine if it makes sense to proceed
		count n = G.numberOfNodes();
		count cn = GandMap.first.numberOfNodes();
		count diff = n - cn;
		repeat = ((diff > 0) &&
				(cn >= MIN_NUM_COMMUNITIES) &&
				((double) diff / (double) n > REL_REPEAT_THRSH)
				); // TODO: last condition: no community becomes too big

		// prepare next iteration if there is one
		if (repeat) {
			G = GandMap.first;
			mapHierarchy.push_back(GandMap.second);
			TRACE("Repeat agglomeration with graph of size " , G.numberOfNodes());
		}
	} while (repeat);

	// vertices of coarsest graph are the clusters
	count cn = G.numberOfNodes();
	Partition zetaCoarse(cn);
	zetaCoarse.allToSingletons();

	// project clustering back to finest graph
	ClusteringProjector projector;
	Partition zeta = projector.projectBackToFinest(zetaCoarse, mapHierarchy,
			graph);

	return zeta;
}

std::string ParallelAgglomerativeClusterer::toString() const {
	std::stringstream strm;
	strm << "ParallelAgglomerativeClusterer";
	return strm.str();
}

} /* namespace NetworKit */
