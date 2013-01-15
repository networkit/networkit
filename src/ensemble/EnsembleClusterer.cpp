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
		// Clustering coreG0(n); // FIXME: project back

		// test if new clustering is better
		qNew = this->qm->getQuality(core, G); // FIXME: quality must be calculated with respect to original graph

		DEBUG("qNew = " << qNew << ", qBest = " << qBest);

		if (qNew > qBest) {
			// new clustering is better
			// FIXME: zetaBest_ = &core;
			qBest = qNew;
			zetaBest_ = &core;

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

Clustering EnsembleClusterer::projectBack(Clustering& zetaCoarse,
		std::vector<NodeMap<node> >& maps, Graph& G0) {

	Clustering zetaFine(G0.numberOfNodes());
	G0.forallNodes([&](node v) {
		node sv = v;
		for (auto map : maps) {
			sv = map[sv];
		}
		cluster sc = zetaCoarse[sv];
		zetaFine.addToCluster(sc, v);
	});

	return zetaFine;
}

Clustering EnsembleClusterer::run2(Graph& G) {

	// sub-algorithms
	ClusterContracter contract;
	RegionGrowingOverlapper overlap;

	// hierarchies
	std::vector<Graph> graph;				// hierarchy of graphs G^{i}
	std::vector<Clustering> clustering;		// hierarchy of core clusterings \zeta^{i}
	std::vector<Clustering> clusteringBack;	// hierarchy of core clusterings projected back to the original graph
	std::vector<NodeMap<node> > map;		// hierarchy of maps M^{i->i+1}
	std::vector<double> quality;			// hierarchy of clustering quality values q^{i} = q(\zeta^{i}, G^{0})

	// other data collections
	std::vector<Clustering>	baseClustering;	// collection of base clusterings, reset in each iteration



	bool repeat;
	int i = -1;	// iteration counter, starts with 0 in loop

	graph.push_back(G); 		// store G^{0}
	Clustering empty(0);
	clusteringBack.push_back(empty); // push a dummy clustering so that clusteringBack[i] contains the clustering projected from G^{i} to G^{0}  (there is no clusteringBack[0])


	do {
		i += 1; 	// increment iteration/hierarchy counter


		INFO("EnsembleClusterer *** ITERATION " << i << " ***");
		baseClustering.clear();

		// *** base clusterers calculate base clusterings ***
		for (auto clusterer : baseClusterers) {
			try {
				baseClustering.push_back(clusterer->run(graph[i]));
				TRACE("created clustering with k=" << baseClustering.back().numberOfClusters());
			} catch (...) {
				ERROR("base clusterer failed with exception.");
				throw std::runtime_error("base clusterer failed.");
			}
		}

		// *** overlap clusters to create core clustering ***
		clustering.push_back(overlap.run(graph[i], baseClustering));

		if (i == 0) {			// first iteration
			// *** calculate quality of first core clustering with respect to first graph ***
			quality.push_back(this->qm->getQuality(clustering[i], graph[i]));
			DEBUG("pushed quality: " << quality.back());

			// *** contract the graph according to core clustering **
			auto con = contract.run(graph[i], clustering[i]);	// returns pair (G^{i+1}, M^{i->i+1})
			graph.push_back(con.first);		// store G^{i+1}
			map.push_back(con.second);		// store M^{i->i+1}

			// new graph created => repeat
			repeat = true;
		} else { 	// other iterations
			clusteringBack.push_back(this->projectBack(clustering[i], map, G));
			quality.push_back(this->qm->getQuality(clusteringBack[i], graph[i]));
			DEBUG("pushed quality: " << quality.back());


			// *** test if new core clustering is better than previous one **
			if (quality[i] > quality[i-1]) {
				DEBUG("quality[" << i << "] = " << quality[i] << " > quality[" << (i-1) << "] = " << quality[i-1]);


				auto con = contract.run(graph[i], clustering[i]);	// returns pair (G^{i+1}, M^{i->i+1})
				graph.push_back(con.first);		// store G^{i+1}
				map.push_back(con.second);		// store M^{i->i+1}

				// new graph created => repeat
				repeat = true;
			} else {
				DEBUG("quality[" << i << "] = " << quality[i] << " <= quality[" << (i-1) << "] = " << quality[i-1]);

				// new core clustering is not better => do not contract according to new core clustering and do not repeat
				repeat = false;
			}
		}

	} while (repeat);

	Clustering zetaCoarse = this->finalClusterer->run(graph[i]); // TODO: check: index correct?
	Clustering zetaFine = this->projectBack(zetaCoarse, map, G);

	return zetaFine;
}

} /* namespace EnsembleClustering */
