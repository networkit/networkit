/*
 * EnsemblePreprocessing.cpp
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "EPP.h"

#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../coarsening/ClusteringProjector.h"
#include "../community/JaccardMeasure.h"
#include "../auxiliary/Log.h"
#include "PLM.h"
#include "PLP.h"

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
	ParallelPartitionCoarsening contracter;
	ClusteringProjector projector;

	// data
	baseClusterings.clear();
	baseClusterings.resize(baseClusterers.size(), Partition(G.upperNodeIdBound())); // collection of base clusterings - fill with empty clustering

	// run base clusterers in parallel
	#pragma omp parallel for
	for (index b = 0; b < baseClusterers.size(); b += 1) {
		// FIXME: initialization of base clusterer?
		baseClusterers.at(b)->run();
		baseClusterings.at(b) = baseClusterers.at(b)->getPartition();
	}

	// create core clustering
	core = this->overlap->run(G, baseClusterings);
	// contract graph according to core clustering
	Graph Gcore;
	std::vector<node> fineToCoarse;
	std::tie(Gcore,fineToCoarse) = contracter.run(G,core);
	// send contracted graph to final clusterer
	// TODO: maybe put this in a private helper function as this could be distracting...
	if (auto tmp = dynamic_cast<PLM*>(this->finalClusterer.get())) {
		DEBUG("final clusterer is PLM");
		this->finalClusterer.reset(new PLM(Gcore, tmp));
	} else if (auto tmp = dynamic_cast<PLP*>(this->finalClusterer.get())) {
		DEBUG("final clusterer is PLP");
		this->finalClusterer.reset(new PLP(Gcore, *tmp));
	}
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


Partition EPP::getCorePartition() const {
	return core;
}

std::vector<Partition> EPP::getBasePartitions() const {
	return baseClusterings;
}

} /* namespace NetworKit */
