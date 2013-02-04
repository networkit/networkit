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

Clustering ClusteringProjector::projectBackToFinest(Clustering& zetaCoarse,
		std::vector<NodeMap<node> >& maps, Graph& Gfinest) {

	std::cout << "# maps: " << maps.size() << std::endl;
	for (auto map : maps) {
		std::cout << "map size: " << map.numberOfNodes() << std::endl;
	}
	std::cout << "zetaCoarse # nodes " << zetaCoarse.numberOfNodes() << std::endl;
	std::cout << "zetaCoarse # clusters " << zetaCoarse.numberOfClusters() << std::endl;


	Clustering zetaFine(Gfinest.numberOfNodes());
	Gfinest.forallNodes([&](node v) {
		node sv = v;
		for (auto map : maps) {
			sv = map.at(sv);
		}
		cluster sc = zetaCoarse[sv];
		zetaFine.addToCluster(sc, v);
	});

	return zetaFine;
}

} /* namespace EnsembleClustering */
