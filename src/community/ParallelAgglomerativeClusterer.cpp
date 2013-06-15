/*
 * ParallelAgglomerativeClusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu),
 *      		Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include "ParallelAgglomerativeClusterer.h"

namespace NetworKit {

ParallelAgglomerativeClusterer::ParallelAgglomerativeClusterer() {
	// TODO Auto-generated constructor stub

}

ParallelAgglomerativeClusterer::~ParallelAgglomerativeClusterer() {
	// TODO Auto-generated destructor stub
}

Clustering ParallelAgglomerativeClusterer::run(Graph& graph) {
	Graph G = graph; // G is the community graph, starts with singletons
	std::vector<NodeMap<node> > mapHierarchy;

	bool repeat = true;
	do {
		// prepare attributes for scoring
		int attrId = G.addEdgeAttribute_double(0.0);

		// perform scoring
		TRACE("before scoring graph of size " << G.numberOfNodes());
		ModularityScoring<double> modScoring(G);
		modScoring.scoreEdges(attrId);

		// FIXME: so far only sequential
		// FIXME: error in second iteration in matching, apparently node
		// degrees are not correctly determined
		// compute matching
		PathGrowingMatcher parMatcher;
		Matching M = parMatcher.run(G);

		// contract graph according to matching, TODO: (and star-like structures)
		MatchingContracter matchingContracter;
		std::pair<Graph, NodeMap<node> > GandMap = matchingContracter.run(G, M);

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
			INFO("Repeat agglomeration with graph of size " << G.numberOfNodes());
		}
	} while (repeat);

	// vertices of coarsest graph are the clusters
	count cn = G.numberOfNodes();
	Clustering zetaCoarse(cn);
	zetaCoarse.allToSingletons();

	// project clustering back to finest graph
	ClusteringProjector projector;
	Clustering zeta = projector.projectBackToFinest(zetaCoarse, mapHierarchy,
			graph);

	return zeta;
}

std::string ParallelAgglomerativeClusterer::toString() const {
	std::stringstream strm;
	strm << "ParallelAgglomerativeClusterer";
	return strm.str();
}

} /* namespace NetworKit */
