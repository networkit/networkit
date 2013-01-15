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

Clustering ClusteringProjector::projectBack(Graph& Gcoarse, Graph& Gfine, NodeMap<node>& fineToCoarse,
		Clustering& zetaCoarse) {

	Clustering zetaFine(Gfine.numberOfNodes());
	// DEBUG
	std::ostringstream oss;	oss << "zeta(" << Gfine.getName() << ")"; zetaFine.setName(oss.str());	//C++??!!
	// DEBUG

	Gfine.forallNodes([&](node v) {
		node sv = fineToCoarse[v];
		cluster cv = zetaCoarse.clusterOf(sv);
		zetaFine[v] = cv;
	});

	return zetaFine;
}

} /* namespace EnsembleClustering */
