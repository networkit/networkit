/*
 * EnsembleClusterer.cpp
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#include "EnsembleClusterer.h"

namespace EnsembleClustering {

EnsembleClusterer::EnsembleClusterer() : Clusterer(),
		baseClusterers(), baseClusterings() {
	this->bestClustering = NULL;
	this->finalClusterer = NULL;
	this->qm = NULL; // TODO: select quality measure
	this->qBest = -1; // TODO: ?
}

EnsembleClusterer::~EnsembleClusterer() {
	// TODO Auto-generated destructor stub
}

bool EnsembleClusterer::isBetterClustering(const Clustering& zeta) {
	// FIXME: double qNew = this->qm->getQuality(zeta);
}

Clustering& EnsembleClusterer::run(Graph& G) {

	// initialize quality measure
	this->qm = new Modularity; // FIXME:initialize in constructor

	ClusteringGenerator clustGen;
	Clustering& zetaSingleton = clustGen.makeSingletonClustering(G);

	this->bestClustering = &zetaSingleton; // TODO: ?

	// store all contracted graphs here
	std::vector<Graph*> contractionHierarchy(log(G.numberOfNodes()), NULL); // TODO: ?

	// base clusterers calculate base clusterings
	for (auto clusterer : baseClusterers) {
		Clustering& zeta = clusterer->run(G);
		baseClusterings.push_back(&zeta);
	}

	// calculate overlap of base clusterings
	RegionGrowingOverlapper overlapper;
	// TODO: Clustering& core = overlapper.run(G, baseClusterings);



	// if new clustering is better than best clustering


	// contract graph according to overlap clustering
	ClusterContracter contracter;
	// TODO: Graph& Gcon = contracter.run(G, core);

	// submit contracted graph to base clusterers again

	// use final clusterer to find clustering of contracted graph
	// TODO: Clustering& zetaFinal = this->finalClusterer->run(GconFinal);

	// project final clustering back to original graph
	// TODO: Clustering& zetaEnsemble = this->projectBack(zetaFinal, GconFinal, G);

	// TODO: implement

	// TODO: return references directly
}

void EnsembleClusterer::addBaseClusterer(Clusterer& baseClusterer) {
	this->baseClusterers.push_back(&baseClusterer);
}

void EnsembleClusterer::setFinalClusterer(Clusterer& clusterer) {
	this->finalClusterer = &clusterer;
}

} /* namespace EnsembleClustering */
