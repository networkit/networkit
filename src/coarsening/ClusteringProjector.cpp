/*
 * ClusteringProjector.cpp
 *
 *  Created on: 07.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringProjector.h"

namespace NetworKit {

ClusteringProjector::ClusteringProjector() {
	// TODO Auto-generated constructor stub

}

ClusteringProjector::~ClusteringProjector() {
	// TODO Auto-generated destructor stub
}

Partition ClusteringProjector::projectBack(Graph& Gcoarse, Graph& Gfine, std::vector<node>& fineToCoarse,
		Partition& zetaCoarse) { // FIXME Partition used like it is supposed to?

	Partition zetaFine(Gfine.numberOfNodes());
	zetaFine.setUpperBound(zetaCoarse.upperBound());
	Gfine.forNodes([&](node v) {
		node sv = fineToCoarse[v];
		index cv = zetaCoarse.subsetOf(sv);
		zetaFine.addToSubset(cv,v);
	});

	return zetaFine;
}

Partition ClusteringProjector::projectBackToFinest(Partition& zetaCoarse,
		std::vector<std::vector<node> >& maps, Graph& Gfinest) {
	if (zetaCoarse.numberOfElements() == Gfinest.numberOfNodes()) {
		return zetaCoarse;
	}

	Partition zetaFine(Gfinest.numberOfNodes());
	zetaFine.setUpperBound(zetaCoarse.upperBound()); // upper bound for ids in zetaFine must be set to upper bound in zetaCoarse, or modularity assertions fail

	// store temporarily coarsest supernode here
	std::vector<node> tempMap(Gfinest.numberOfNodes());
	// TODO: add NodeMap.setAll, NodeMap.identity...
	// initialize to identity
	Gfinest.parallelForNodes([&](node v){
		tempMap[v] = v;
	});

	// find coarsest supernode for each node
	for (auto iter = maps.begin(); iter != maps.end(); ++iter) {
		Gfinest.parallelForNodes([&](node v){
			tempMap[v] = (* iter)[tempMap[v]];
		});
	}


	// set clusters for fine nodes
	Gfinest.parallelForNodes([&](node v) {
		index sc = zetaCoarse[tempMap[v]];
		zetaFine.addToSubset(sc,v);//zetaFine[v] = sc;
	});

	return zetaFine;
}

Partition ClusteringProjector::projectCoarseGraphToFinestClustering(Graph& Gcoarse, Graph& Gfinest, std::vector<std::vector<node> >& maps) {

	Partition zeta(Gfinest.numberOfNodes());

	// store temporarily coarsest supernode here
	std::vector<node> super(Gfinest.numberOfNodes());
	// initialize to identity
	Gfinest.parallelForNodes([&](node v){
		super[v] = v;
	});


	// find coarsest supernode for each node
	for (auto iter = maps.begin(); iter != maps.end(); ++iter) {
		Gfinest.parallelForNodes([&](node v){
			super[v] = (* iter)[super[v]];
		});
	}

	// assign super node id as cluster id
	Gfinest.parallelForNodes([&](node v) {
		//zeta.addToSubset(super[v],v);
		zeta[v] = super[v];
	});

	// Partition not used properly...set upperBound remotely...	
	index maxID = 0;
	zeta.forEntries([&](index e, index s) {
		if (s > maxID) {
			maxID = s;
		}
	});
	zeta.setUpperBound(maxID);

	DEBUG("number of clusters in projected clustering: " , zeta.numberOfSubsets());
	DEBUG("number of nodes in coarse graph: " , Gcoarse.numberOfNodes());
	assert (zeta.numberOfSubsets() == Gcoarse.numberOfNodes());

	return zeta;

}
/* FIXME: to be deleted?
Partition ClusteringProjector::projectBack(Graph& Gcoarse, Graph& Gfine, std::vector<node>& fineToCoarse, Partition& zetaCoarse) {

	Partition zetaFine(Gfine.numberOfNodes());

	Gfine.forNodes([&](node v) {
		node sv = fineToCoarse[v];
		index cv = zetaCoarse[sv];//zetaCoarse.clusterOf(sv);
		zetaFine.addToSubset(cv,v); //zetaFine[v] = cv;
	});

	return zetaFine;
}*/

} /* namespace NetworKit */


