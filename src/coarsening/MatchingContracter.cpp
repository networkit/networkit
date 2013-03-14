/*
 * MatchingContracter.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "MatchingContracter.h"

namespace EnsembleClustering {

MatchingContracter::MatchingContracter() {
	// TODO Auto-generated constructor stub

}

MatchingContracter::~MatchingContracter() {
	// TODO Auto-generated destructor stub
}

std::pair<Graph, NodeMap<node> > MatchingContracter::run(Graph& G, Matching& M) {
	count n = G.numberOfNodes();
	count cn = n - M.matchingSize();
	Graph cG(cn);

	// compute map: old ID -> new coarse ID
	index idx = 0;
	NodeMap<node> mapFineToCoarse(n);
	G.forNodes([&](node v) { // TODO: difficult in parallel
		index mate = M.mate(v);
		if ((mate == none) || (v < mate)) {
			// vertex is carried over to the new level
			mapFineToCoarse[v] = idx;
			++idx;
		}
		else {
			// vertex is not carried over, receives ID of mate
			mapFineToCoarse[v] = mapFineToCoarse[mate];
		}
	});

//	for (node v = 0; v < n; ++v) {
//		std::cout << v << " maps to " << mapFineToCoarse[v] << std::endl;
//	}
//	std::cout << "matching size: " << M.matchingSize() << std::endl;

	G.forNodes([&](node v) { // TODO: difficult in parallel
		G.forNeighborsOf(v, [&](node u) {
			node cv = mapFineToCoarse[v];
			node cu = mapFineToCoarse[u];
			edgeweight ew = G.weight(v, u);
			cG.setWeight(cv, cu, cG.weight(cv, cu) + ew);
		});
	});

	return std::make_pair(cG, mapFineToCoarse);
}

} /* namespace EnsembleClustering */
