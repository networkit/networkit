/*
 * EnsembleClusterer.cpp
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#include "EnsembleClusterer.h"

namespace EnsembleClustering {

EnsembleClusterer::EnsembleClusterer() : Clusterer() {
	this->finalClusterer = NULL;
	this->qm = NULL;
}

EnsembleClusterer::~EnsembleClusterer() {
	// TODO Auto-generated destructor stub
}

void EnsembleClusterer::setQualityMeasure(QualityMeasure& qm) {
	this->qm = &qm;
}


void EnsembleClusterer::addBaseClusterer(Clusterer& base) {
	this->baseClusterers.push_back(&base);
}

void EnsembleClusterer::setFinalClusterer(Clusterer& final) {
	this->finalClusterer = &final;
}

Clustering EnsembleClusterer::run(Graph& G) {

	// sub-algorithms
	ClusterContracter contracter;
	RegionGrowingOverlapper overlapper;
	// data
	Clustering* zetaBest = NULL;
	std::vector<Clustering> baseClusterings;
	double qNew;	// quality of new clustering
	double qBest;	// quality of best clustering
	bool repeat; 	// loop condition
	Graph* Gbest = NULL;	// contracted graph belonging to the best clustering

	// store all clusterings/contracted graphs here
	double logn = log(G.numberOfNodes());
	std::vector<Graph> graphHierarchy(logn, NULL);
	std::vector<Clustering> clusteringHierarchy(logn, NULL);

	qBest = -1; // => repeat loop at least once
	do {
		baseClusterings.clear();

		// base clusterers calculate base clusterings
		for (auto clusterer : baseClusterers) {
			Clustering zeta = clusterer->run(G);
			baseClusterings.push_back(zeta);
		}

		// calculate overlap of base clusterings
		Clustering core = overlapper.run(G, baseClusterings);
		// store clustering & graph
		clusteringHierarchy.push_back(core);
		graphHierarchy.push_back(G);

		if (zetaBest == NULL) {
			// in first iteration
			zetaBest = &core;
			qBest = -1; // contract graph at least once
			repeat = true; // repeat at least once
		} else {
			// after first iteration
			// test if new clustering is better
			qNew = this->qm->getQuality(core, G); // FIXME: quality must be calculated with respect to original graph - project back
			DEBUG("qNew = " << qNew << ", qBest = " << qBest);
			if (qNew > qBest) {
				// new clustering is better
				zetaBest = &core;
				qBest = qNew;
				// contract graph according to overlap clustering
				GraphContraction con = contracter.run(G, core);
				Graph Gcon = con.getCoarseGraph();
				// work on contracted graph
				Gbest = &Gcon;
				repeat = true; // submit contracted graph to base clusterers again
			} else {
				// new clustering is not better
				repeat = false;
			}
		}
		// repeat while new clustering is better than old clustering
	} while (qNew > qBest);

	// work on contracted graph for best clustering

	// use final clusterer to find clustering of contracted graph
	Clustering zetaFinal = this->finalClusterer->run(*Gbest); // FIXME: Gbest can be NULL


	// project final clustering back to original graph
	// TODO: Clustering& zetaEnsemble = this->projectBack(zetaFinal, GconFinal, G);


	return zetaFinal;
}


} /* namespace EnsembleClustering */
