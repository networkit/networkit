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

	// data storage
	std::vector<Clustering> clusteringHierarchy;
	std::vector<GraphContraction> contractionHierarchy;

	std::vector<Clustering> baseClusterings;

	// temporaries
	Graph* G_ = &G;				// the (fine) graph which is to be clustered by the base clusterers
	Graph* Gcoarse_ = NULL;
	Graph* Gbest_ = NULL;		// contracted graph belonging to the best clustering - TODO: needed?
	Clustering* zetaBest_ = NULL;

	// values
	double qNew	= -1;		// quality of new clustering
	double qBest = -1;		// quality of best clustering => repeat loop at least once
	int nIter = 0;			// number of iterations
	bool repeat;	// loop condition


	do {
		nIter += 1;
		INFO("EnsembleClusterer *** ITERATION " << nIter << " ***");
		baseClusterings.clear();

		// base clusterers calculate base clusterings
		for (auto clusterer : baseClusterers) {
			try {
				Clustering zeta = clusterer->run(G);
				baseClusterings.push_back(zeta);
				TRACE("created clustering with k=" << zeta.numberOfClusters());
			} catch (...) {
				ERROR("base clusterer failed with exception.");
				throw std::runtime_error("base clusterer failed.");
			}
		}

		// calculate overlap of base clusterings
		Clustering core = overlapper.run(G, baseClusterings);
		clusteringHierarchy.push_back(core);

		DEBUG("created core clustering with k=" << core.numberOfClusters());

		// test if new clustering is better
		Clustering coreG0
		qNew = this->qm->getQuality(core, G); // FIXME: quality must be calculated with respect to original graph - project back

		DEBUG("qNew = " << qNew << ", qBest = " << qBest);

		if (qNew > qBest) {
			// new clustering is better
			zetaBest_ = &core;
			qBest = qNew;
			// contract graph according to overlap clustering
			contractionHierarchy.push_back(contracter.run(G, core)); 			// store contraction in hierarchy

			Gcoarse_ = &(contractionHierarchy.back().getCoarseGraph());

			// use this
//			Graph* Gcon = &(contractionHierarchy[i].getCoarseGraph());
//
//			Graph Gcon = con.getCoarseGraph();
//
//			Graph* Gcon = new Graph(con.getCoarseGraph());
//
//			Graph* Gcon = &con.getCoarseGraph();

			// DEBUG
			std::ostringstream oss;	oss << "G^" << nIter; Gcoarse_->setName(oss.str());	//C++??!!
			// DEBUG
			INFO("contracted graph created: " << Gcoarse_->toString());
			// work on contracted graph
			Gbest_ = Gcoarse_;
			// FIXME: G = Gcon;	// TODO: correct?
			// FIXME: DEBUG("G is now " << G.toString());
			repeat = true;
		} else {
			// new clustering is not better
			repeat = false;
		}
		// repeat while new clustering is better than old clustering
		// FIXME: set new

	} while (repeat);

	DEBUG("LOOP EXIT");

	// work on contracted graph for best clustering

	// use final clusterer to find clustering of contracted graph
	assert (Gbest_ != NULL);
	DEBUG("Gbest_: " << Gbest_->toString());
	Clustering zetaFinal = this->finalClusterer->run(*Gbest_); // FIXME: LabelPropagation does not terminate


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
