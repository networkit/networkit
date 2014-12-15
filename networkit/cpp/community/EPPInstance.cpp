/*
 * EnsemblePreprocessing.cpp
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "EPPInstance.h"

#include "../coarsening/ClusterContractor.h"
#include "../coarsening/ClusteringProjector.h"
#include "../community/JaccardMeasure.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

EPPInstance::EPPInstance(const Graph& G, count ensembleSize) : CommunityDetectionAlgorithm(G), ensembleSize(ensembleSize) {
}


void EPPInstance::run() {
	INFO("STARTING EnsemblePreprocessing on G=" , G.toString());

	// fixed sub-algorithms
	std::vector<PLP> baseClusterers; //!< ensemble of base clusterers
	HashingOverlapper overlap; //!< clustering overlap algorithm


	ClusterContractor contracter;
	ClusteringProjector projector;

	// data
	INFO("init base clusterings");
	baseClusterings.clear();
	baseClusterings.resize(ensembleSize, Partition(G.upperNodeIdBound())); // collection of base clusterings - fill with empty clustering

	//
	INFO("init base PLPs");
	for (index b = 0; b < ensembleSize; b += 1) {
		baseClusterers.push_back(PLP(G));
	}

	INFO("run base PLPs");
	// run base clusterers in parallel
	#pragma omp parallel for
	for (index b = 0; b < ensembleSize; b += 1) {
		baseClusterers.at(b).run();
		baseClusterings.at(b) = baseClusterers.at(b).getPartition();
	}

	// ANALYSIS
	if (CALC_DISSIMILARITY) {
		JaccardMeasure dm;
		double dissimilaritySum = 0.0;
		for (index b = 0; b < ensembleSize; b += 1) {
			for (index c = b + 1; c < ensembleSize; c += 1) {
				double d = dm.getDissimilarity(G, baseClusterings.at(b), baseClusterings.at(c));
				dissimilaritySum += d;
			}
		}
		double avgDissimilarity = dissimilaritySum / (baseClusterings.size() * (baseClusterings.size() - 1) / 2.0);
		std::cout << "[INFO] avg. base clustering dissimilarity: " << avgDissimilarity << std::endl;
	}
	//

	INFO("overlap");
	// create core clustering
	core = overlap.run(G, baseClusterings);
	// contract graph according to core clustering
	INFO("coarsening");
	std::pair<Graph, std::vector<node> > contraction = contracter.run(G, core);
	Graph Gcore = contraction.first;
	std::vector<node> fineToCoarse = contraction.second;
	// send contracted graph to final clusterer
	PLM finalClusterer(Gcore, true);	//!< final clustering algorithm: PLMR
	finalClusterer.run();
	Partition finalCoarse = finalClusterer.getPartition();

	// project clustering of contracted graph back to original graph
	Partition final = projector.projectBack(Gcore, G, fineToCoarse, finalCoarse);
	// return clustering
	result = std::move(final);
	hasRun = true;
}

std::string EPPInstance::toString() const {
	std::stringstream strm;
	strm << "EPPInstance(" << ensembleSize << ")";
	return strm.str();
}


Partition EPPInstance::getCorePartition() const {
	return core;
}

std::vector<Partition> EPPInstance::getBasePartitions() const {
	return baseClusterings;
}

} /* namespace NetworKit */
