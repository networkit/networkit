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
	// std::vector<GraphContraction> contractionHierarchy;
	std::vector<Graph> graphHierarchy;
	std::vector<NodeMap<node> > mapHierarchy;

	std::vector<Clustering> baseClusterings;

	// temporaries
	Graph* G0_ = &G;			// points to the input graph G^0
	Graph* G_ = &G;				// points to the current (fine) graph which is to be clustered by the base clusterers
	Graph* Gcoarse_ = NULL;		// points to the
	// Graph* Gbest_ = NULL;		// contracted graph belonging to the best clustering - TODO: needed?
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
				Clustering zeta = clusterer->run(*G_);
				baseClusterings.push_back(zeta);
				TRACE("created clustering with k=" << zeta.numberOfClusters());
			} catch (...) {
				ERROR("base clusterer failed with exception.");
				throw std::runtime_error("base clusterer failed.");
			}
		}

		// calculate overlap of base clusterings
		Clustering core = overlapper.run(*G_, baseClusterings);
		clusteringHierarchy.push_back(core);

		DEBUG("created core clustering with k=" << core.numberOfClusters());

		// TODO: project core clustering back to G
		Clustering coreG0; // FIXME: project back

		// test if new clustering is better
		qNew = this->qm->getQuality(coreG0, G); // FIXME: quality must be calculated with respect to original graph

		DEBUG("qNew = " << qNew << ", qBest = " << qBest);

		if (qNew > qBest) {
			// new clustering is better
			// FIXME: zetaBest_ = &core;
			qBest = qNew;
			zetaBest_ = &coreG0;

			std::pair<Graph, NodeMap<node> > contraction = contracter.run(G, core);
			// contract graph according to overlap clustering
			graphHierarchy.push_back(contraction.first); 	// store coarse graph in hierarchy
			mapHierarchy.push_back(contraction.second);		// store fine->coarse node map in hierarchy

			Gcoarse_ = &(graphHierarchy.back());

			// DEBUG
			std::ostringstream oss;	oss << "G^" << nIter; Gcoarse_->setName(oss.str());	//C++??!!
			// DEBUG
			INFO("contracted graph created: " << Gcoarse_->toString());

			// work on coarse graph now
			G_ = Gcoarse_;
			// if new coarse graph has been created, repeat
			repeat = true;

		} else {
			// new clustering is not better
			repeat = false;
		}
		// repeat while new clustering is better than old clustering
	} while (repeat);

	DEBUG("LOOP EXIT");

	// work on contracted graph for best clustering

	// use final clusterer to find clustering of contracted graph
	assert (Gbest_ != NULL);
	DEBUG("Gbest_: " << Gbest_->toString());
	Clustering zetaFinal = this->finalClusterer->run(*Gbest_); // FIXME: LabelPropagation does not terminate
	zetaFinal.setName("zetaFinal");


	// project final clustering back to original graph
	ClusteringProjector projector;
	int h = graphHierarchy.size() - 1; // coarsest contraction in the hierarchy

	// DEBUG
	DEBUG("graph hierarchy size: " << h);
	assert (h > 0);
	// DEBUG

	Clustering zetaBack = zetaFinal;
	for (int i = h; i >= 1; --i) {
		zetaBack = projector.projectBack(graphHierarchy.at(i), graphHierarchy.at(i - 1), mapHierarchy.at(i - 1), zetaBack);
	}

	DEBUG("returning zetaBack: " << zetaBack.getName());
	return zetaBack;
}


} /* namespace EnsembleClustering */
