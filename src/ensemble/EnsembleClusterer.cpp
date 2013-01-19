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
	// DEBUG
	INFO("STARTING EnsembleClusterer on G=" << G.toString());
	// DEBUG

	// sub-algorithms
	ClusterContracter contract;
	RegionGrowingOverlapper overlap;
	ClusteringProjector project;

	// hierarchies
	std::vector<Graph> graph;				// hierarchy of graphs G^{i}
	std::vector<Clustering> clustering;		// hierarchy of core clusterings \zeta^{i}
	std::vector<Clustering> clusteringBack;	// hierarchy of core clusterings projected back to the original graph
	std::vector<NodeMap<node> > map;		// hierarchy of maps M^{i->i+1}
	std::vector<double> quality;			// hierarchy of clustering quality values q^{i} = q(\zeta^{i}, G^{0})

	// other data collections
	std::vector<Clustering>	baseClustering(baseClusterers.size(), Clustering(0));	// collection of base clusterings - fill with empty clustering


	// DEBUG
	GraphIO graphio;
	// DEBUG


	bool repeat;
	int i = -1;	// iteration counter, starts with 0 in loop

	graph.push_back(G); 		// store G^{0}
	Clustering empty(0);
	clusteringBack.push_back(empty); // push a dummy clustering so that clusteringBack[i] contains the clustering projected from G^{i} to G^{0}  (there is no clusteringBack[0])


	do {
		i += 1; 	// increment iteration/hierarchy counter


		INFO("EnsembleClusterer *** ITERATION " << i << " ***");

		// *** base clusterers calculate base clusterings ***
		#pragma omp parallel for
		for (int b = 0; b < baseClusterers.size(); b += 1) {	// TODO: run base clusterers in parallel
			try {
				baseClustering[b] = baseClusterers[b]->run(graph[i]);	// TODO: is push_back a critical section?
				// DEBUG
				DEBUG("created base clustering: k=" << baseClustering[b].numberOfClusters());
				if (baseClustering[b].isOneClustering(graph[i])) {
					WARN("base clusterer created 1-clustering");
				}
				// DEBUG
			} catch (...) {
				ERROR("base clusterer failed with exception.");
				throw std::runtime_error("base clusterer failed.");
			}
		}

		// *** overlap clusters to create core clustering ***
		clustering.push_back(overlap.run(graph[i], baseClustering));
		// DEBUG
		DEBUG("created core clustering: k=" << clustering[i].numberOfClusters());

		// DEBUG

		if (i == 0) {			// first iteration
			// *** calculate quality of first core clustering with respect to first graph ***
			quality.push_back(this->qm->getQuality(clustering[i], graph[i]));
			DEBUG("pushed quality: " << quality.back());

			// *** contract the graph according to core clustering **
			auto con = contract.run(graph[i], clustering[i]);	// returns pair (G^{i+1}, M^{i->i+1})
			graph.push_back(con.first);		// store G^{i+1}
			map.push_back(con.second);		// store M^{i->i+1}

			//DEBUG
			DEBUG("contracted graph G^" << (i+1) << " created: " << graph.back().toString());
			//DEBUG

			// new graph created => repeat
			repeat = true;
		} else { 	// other iterations
			clusteringBack.push_back(project.projectBackToFinest(clustering[i], map, G));
			assert (clustering[i].numberOfClusters() == clusteringBack[i].numberOfClusters());
			// DEBUG
			DEBUG("created projected clustering: k=" << clusteringBack[i].numberOfClusters());
			// DEBUG

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
	Clustering zetaFine = project.projectBackToFinest(zetaCoarse, map, G);

	return zetaFine;
}

} /* namespace EnsembleClustering */
