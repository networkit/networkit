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
	// FIXME: LabelPropagation does not terminate on contracted graph
	DEBUG("starting EnsembleClusterer on graph G with n=" << G.numberOfNodes() << " m=" << G.numberOfEdges());



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
	std::vector<Graph> graphHierarchy;
	std::vector<Clustering> clusteringHierarchy;
	std::vector<GraphContraction> contractionHierarchy; // TODO: store

	int nIter = 0;
	qBest = -1; // => repeat loop at least once
	do {
		nIter += 1;
		INFO("\n EnsembleClusterer *** ITERATION " << nIter << " ***");
		baseClusterings.clear();

		// store current graph to be clustered by ensemble
		graphHierarchy.push_back(G);

		// base clusterers calculate base clusterings
		for (auto clusterer : baseClusterers) {
			try {
				Clustering zeta = clusterer->run(G);
				baseClusterings.push_back(zeta);
				DEBUG("created clustering with k=" << zeta.numberOfClusters());
			} catch (...) {
				ERROR("base clusterer failed.");
				throw std::runtime_error("base clusterer failed.");
			}
		}

		// calculate overlap of base clusterings
		Clustering core = overlapper.run(G, baseClusterings);
		DEBUG("created core clustering with k=" << core.numberOfClusters());

		// store core clustering
		clusteringHierarchy.push_back(core);

		// FIXME: contract graph also in first iteration

		if (zetaBest == NULL) {
			// in first iteration
			zetaBest = &core;
			qBest = -1; // contract graph at least once
			repeat = true; // repeat at least once
		}

		// test if new clustering is better
		qNew = this->qm->getQuality(core, G); // FIXME: quality must be calculated with respect to original graph - project back
		DEBUG("qNew = " << qNew << ", qBest = " << qBest);
		if (qNew > qBest) {
			// new clustering is better
			zetaBest = &core;
			qBest = qNew;
			// contract graph according to overlap clustering
			GraphContraction con = contracter.run(G, core);
			contractionHierarchy.push_back(con); 			// store contraction in hierarchy

			Graph Gcon = con.getCoarseGraph();
			DEBUG("contracted graph: n=" << Gcon.numberOfNodes() << " m=" << Gcon.numberOfEdges());
			// work on contracted graph
			Gbest = &Gcon;
			repeat = true; // submit contracted graph to base clusterers again
		} else {
			// new clustering is not better
			repeat = false;
		}
		// repeat while new clustering is better than old clustering
	} while (qNew > qBest);

	// work on contracted graph for best clustering

	// use final clusterer to find clustering of contracted graph
	assert (Gbest != NULL);
	DEBUG("Gbest graph: n=" << Gbest->numberOfNodes() << " m=" << Gbest->numberOfEdges());
	Clustering zetaFinal = this->finalClusterer->run(*Gbest); // FIXME: LabelPropagation does not terminate


	// project final clustering back to original graph
	ClusteringProjector projector;
	int h = contractionHierarchy.size() - 1; // coarsest contraction in the hierarchy
	Clustering zeta = zetaFinal;
	for (int i = h; i >= 0; --i) {
		zeta = projector.projectBack(contractionHierarchy.at(i), zeta); // TODO: check
	}

	// FIXME: clustering returned does not belong to G
	return zeta;
}


} /* namespace EnsembleClustering */
