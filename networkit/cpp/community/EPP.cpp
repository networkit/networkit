/*
 * EnsemblePreprocessing.cpp
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "EPP.h"

#include "../coarsening/ClusterContractor.h"
#include "../coarsening/ClusteringProjector.h"
#include "../community/JaccardMeasure.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

EPP::EPP(const Graph& G) : CommunityDetectionAlgorithm(G) {
	this->finalClusterer = NULL;
	this->overlap = NULL;
}


void EPP::addBaseClusterer(std::unique_ptr<CommunityDetectionAlgorithm>& base) {
	this->baseClusterers.push_back(std::move(base));
}

void EPP::setFinalClusterer(std::unique_ptr<CommunityDetectionAlgorithm>& final) {
	this->finalClusterer = std::move(final);
}

void EPP::setOverlapper(std::unique_ptr<Overlapper>& overlap) {
	this->overlap = std::move(overlap);
}

void EPP::run() {
	INFO("STARTING EnsemblePreprocessing on G=" , G.toString());

	// fixed sub-algorithms
	ClusterContractor contracter;
	ClusteringProjector projector;

	// data
	std::vector<Partition> baseClusterings(baseClusterers.size(), Partition(0)); // collection of base clusterings - fill with empty clustering

	// run base clusterers in parallel
	#pragma omp parallel for
	for (index b = 0; b < baseClusterers.size(); b += 1) {
		// FIXME: initialization of base clusterer?
		baseClusterers.at(b)->run();
		baseClusterings.at(b) = baseClusterers.at(b)->getPartition();
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
	Partition core = this->overlap->run(G, baseClusterings);
	// contract graph according to core clustering
	std::pair<Graph, std::vector<node> > contraction = contracter.run(G, core);
	Graph Gcore = contraction.first;
	std::vector<node> fineToCoarse = contraction.second;
	// send contracted graph to final clusterer
	// FIXME: initialization of final clusterer?
	this->finalClusterer->run();
	Partition finalCoarse = this->finalClusterer->getPartition();

	// project clustering of contracted graph back to original graph
	Partition final = projector.projectBack(Gcore, G, fineToCoarse, finalCoarse);
	// return clustering
	result = std::move(final);
	hasRun = true;
}

std::string EPP::toString() const {
	std::stringstream strm;
	strm << "EnsemblePreprocessing(" << "base=" << this->baseClusterers.front()->toString() << ",ensemble=" << this->baseClusterers.size() << ",final=" << this->finalClusterer->toString() << ")";
	return strm.str();
}


} /* namespace NetworKit */
