/*
 * ClusteringProjector.cpp
 *
 *  Created on: 07.01.2013
 *      Author: cls
 */

#include "ClusteringProjector.h"

namespace EnsembleClustering {

ClusteringProjector::ClusteringProjector() {
	// TODO Auto-generated constructor stub

}

ClusteringProjector::~ClusteringProjector() {
	// TODO Auto-generated destructor stub
}

Clustering ClusteringProjector::projectBack(GraphContraction& contraction,
		Clustering& zetaCoarse) {

	Graph Gfine = contraction.getFineGraph();
	Graph Gcoarse = contraction.getCoarseGraph();
	auto fineToCoarse = contraction.getFineToCoarseMap();

	Clustering zetaFine(Gfine.numberOfNodes());

	Gfine.forallNodes([&](node v) {
		node sv = fineToCoarse[v];
		cluster cv = zetaCoarse.clusterOf(sv);
		zetaFine[v] = cv;
	});

	return zetaFine;
}

} /* namespace EnsembleClustering */
