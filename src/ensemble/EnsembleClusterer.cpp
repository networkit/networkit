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
	G.setName("G^0");
	INFO("starting EnsembleClusterer on graph G: " << G.toString());



	// sub-algorithms
	ClusterContracter contracter;
	RegionGrowingOverlapper overlapper;
	// data
	Clustering* zetaBest = NULL;
	std::vector<Clustering> baseClusterings;
	double qNew;	// quality of new clustering
	double qBest;	// quality of best clustering
	Graph* Gbest = NULL;	// contracted graph belonging to the best clustering

	// store all clusterings/contracted graphs here
	std::vector<Graph> graphHierarchy;
	std::vector<Clustering> clusteringHierarchy;
	std::vector<GraphContraction> contractionHierarchy; // TODO: store

	int nIter = 0;
	qBest = -1; // => repeat loop at least once
	bool repeat;	// loop condition
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
		INFO("created core clustering with k=" << core.numberOfClusters());

		// store core clustering
		clusteringHierarchy.push_back(core);


		if (zetaBest == NULL) {
			// in first iteration
			zetaBest = &core;
			qBest = -1; // contract graph at least once
			repeat = true; // repeat at least once
		}

		// test if new clustering is better
		qNew = this->qm->getQuality(core, G); // FIXME: quality must be calculated with respect to original graph - project back

		INFO("qNew = " << qNew << ", qBest = " << qBest);

		if (qNew > qBest) {
			// new clustering is better
			zetaBest = &core;
			qBest = qNew;
			// contract graph according to overlap clustering
			GraphContraction con = contracter.run(G, core);
			contractionHierarchy.push_back(con); 			// store contraction in hierarchy


			Graph Gcon = con.getCoarseGraph();
			// DEBUG
			std::ostringstream oss;	oss << "G^" << nIter; Gcon.setName(oss.str());	//C++??!!
			// DEBUG
			INFO("contracted graph created: " << Gcon.toString());
			// work on contracted graph
			Gbest = &Gcon;
			G = Gcon;	// TODO: correct?
			DEBUG("G is now" << G.toString());
			repeat = true;
		} else {
			// new clustering is not better
			repeat = false;
		}
		// repeat while new clustering is better than old clustering
		// FIXME: set new

	} while (repeat);

	// work on contracted graph for best clustering

	// use final clusterer to find clustering of contracted graph
	assert (Gbest != NULL);
	INFO("Gbest graph: n=" << Gbest->numberOfNodes() << " m=" << Gbest->numberOfEdges());
	Clustering zetaFinal = this->finalClusterer->run(*Gbest); // FIXME: LabelPropagation does not terminate


	// project final clustering back to original graph
	ClusteringProjector projector;
	int h = contractionHierarchy.size() - 1; // coarsest contraction in the hierarchy
	Clustering zeta = zetaFinal;
	for (int i = h; i >= 0; --i) {
		zeta = projector.projectBack(contractionHierarchy.at(i), zeta);
	}

	return zeta;
}


} /* namespace EnsembleClustering */
