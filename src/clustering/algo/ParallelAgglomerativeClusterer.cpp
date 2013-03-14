/*
 * ParallelAgglomerativeClusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
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
	Graph G = graph; // G is the community graph, starts with singletons
	std::vector<NodeMap<node> > mapHierarchy;

	bool repeat = true;
	do {
		// prepare attributes for scoring
		int attrId = G.addEdgeAttribute_double(0.0);

		// perform scoring
		ModularityScoring<double> modScoring(G);
		modScoring.scoreEdges(attrId);

		// compute matching
		ParallelMatcher parMatcher(attrId);
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

} /* namespace EnsembleClustering */
