/*
 * EnsemblePreprocessing.cpp
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "EPP.h"


#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"
#include "../clustering/JaccardMeasure.h"

namespace NetworKit {

EPP::EPP() : Clusterer() {
	this->finalClusterer = NULL;
	this->overlap = NULL;
}


void EPP::addBaseClusterer(Clusterer& base) {
	this->baseClusterers.push_back(&base);
}

void EPP::setFinalClusterer(Clusterer& final) {
	this->finalClusterer = &final;
}

void EPP::setOverlapper(Overlapper& overlap) {
	this->overlap = &overlap;
}

Clustering EPP::run(Graph& G) {
	INFO("STARTING EnsemblePreprocessing on G=" << G.toString());

	// fixed sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// data
	std::vector<Clustering> baseClusterings(baseClusterers.size(), Clustering(0)); // collection of base clusterings - fill with empty clustering

	// run base clusterers in parallel
	#pragma omp parallel for
	for (index b = 0; b < baseClusterers.size(); b += 1) {
		baseClusterings.at(b) = baseClusterers.at(b)->run(G);
	}

	// ANALYSIS
	if (CALC_DISSIMILARITY) {
		JaccardMeasure dm;
		double dissimilaritySum = 0.0;
		for (index b = 0; b < baseClusterings.size(); b += 1) {
			for (index c = b + 1; c < baseClusterings.size(); c += 1) {
				double d = dm.getDissimilarity(G, baseClusterings.at(b), baseClusterings.at(c));
				dissimilaritySum += d;
			}
		}
		double avgDissimilarity = dissimilaritySum / (baseClusterings.size() * (baseClusterings.size() - 1) / 2.0);
		std::cout << "[INFO] avg. base clustering dissimilarity: " << avgDissimilarity << std::endl;
	}
	//

	// create core clustering
	Clustering core = this->overlap->run(G, baseClusterings);
	// contract graph according to core clustering
	std::pair<Graph, NodeMap<node> > contraction = contracter.run(G, core);
	Graph Gcore = contraction.first;
	NodeMap<node> fineToCoarse = contraction.second;
	// send contracted graph to final clusterer
	Clustering finalCoarse = this->finalClusterer->run(Gcore);

	// project clustering of contracted graph back to original graph
	Clustering final = projector.projectBack(Gcore, G, fineToCoarse, finalCoarse);
	// return clustering
	return final;
}

std::string EPP::toString() const {
	std::stringstream strm;
	strm << "EnsemblePreprocessing(" << "base=" << this->baseClusterers.front()->toString() << ",ensemble=" << this->baseClusterers.size() << ",final=" << this->finalClusterer->toString() << ")";
	return strm.str();
}

} /* namespace NetworKit */
