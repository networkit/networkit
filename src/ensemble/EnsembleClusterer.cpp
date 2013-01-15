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

Clustering EnsembleClusterer::run(Graph& G) {

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
				// DEBUG
				DEBUG("created clustering with k=" << baseClustering.back().numberOfClusters());
				if (baseClustering.back().isOneClustering(graph[i])) {
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

			if (graph[i+1].numberOfEdges() == 0) {
				double vw = 0.0;
				G.forallNodes([&](node v) {
					vw += graph[i+1].weight(v);
				});
				DEBUG("graph has no edges, total node (self-loop) weight is:" << vw);
			}
			//DEBUG

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
