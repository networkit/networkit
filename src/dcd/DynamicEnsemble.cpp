/*
 * DynamicEnsemble.cpp
 *
 *  Created on: 19.06.2013
 *      Author: cls
 */

#include "DynamicEnsemble.h"
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"

namespace NetworKit {

DynamicEnsemble::DynamicEnsemble() : overlapAlgo(NULL), finalAlgo(NULL) {
	// TODO Auto-generated constructor stub

}

DynamicEnsemble::~DynamicEnsemble() {
	// TODO Auto-generated destructor stub
}

void DynamicEnsemble::setGraph(Graph& G) {
	if (!G.isEmpty()) {
		throw std::runtime_error("G is not an empty graph. Currently, it is assumed that this algorithm is initialized with an empty graph, which is then constructed incrementally");
	}

	this->G = &G;
	// also set the graph for all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->setGraph(G);
	}
}

Clustering DynamicEnsemble::run() {
	INFO("STARTING DynamicEnsemble on G=" << G->toString());
	this->runtime.start(); // start timer

	// fixed sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// data
	std::vector<Clustering> baseClusterings(baseAlgos.size(), Clustering(0)); // collection of base clusterings - fill with empty clustering

	// run base clusterers in parallel
	#pragma omp parallel for
	for (index b = 0; b < baseAlgos.size(); b += 1) {
		baseClusterings.at(b) = baseAlgos.at(b)->run();
	}

	// create core clustering
	Clustering core = this->overlapAlgo->run(*G, baseClusterings);
	// contract graph according to core clustering
	std::pair<Graph, NodeMap<node> > contraction = contracter.run(*G, core);
	Graph Gcore = contraction.first;
	NodeMap<node> fineToCoarse = contraction.second;
	// send contracted graph to final clusterer
	Clustering finalCoarse = this->finalAlgo->run(Gcore);

	// project clustering of contracted graph back to original graph
	Clustering final = projector.projectBack(Gcore, *G, fineToCoarse, finalCoarse);

	this->runtime.stop(); // stop timer
	this->timerHistory.push_back(this->runtime.elapsedMilliseconds());
	// return clustering
	return final;
}

std::string DynamicEnsemble::toString() const {
	return "DynamicEnsemble";
}

void DynamicEnsemble::onNodeAddition(node u) {
	// delegate event to all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->onNodeAddition(u);
	}
}

void DynamicEnsemble::onNodeRemoval(node u) {
	// delegate event to all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->onNodeRemoval(u);
	}
}

void DynamicEnsemble::onEdgeAddition(node u, node v, edgeweight w) {
	// delegate event to all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->onEdgeAddition(u, v, w);
	}
}

void DynamicEnsemble::onEdgeRemoval(node u, node v, edgeweight w) {
	// delegate event to all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->onEdgeRemoval(u, v, w);
	}
}

void DynamicEnsemble::onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) {
	// delegate event to all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->onWeightUpdate(u, v, wOld, wNew);
	}
}

void DynamicEnsemble::onTimeStep() {
	// delegate event to all base algorithms
	for (DynamicCommunityDetector* algo : baseAlgos) {
		algo->onTimeStep();
	}
}

void DynamicEnsemble::addBaseAlgorithm(DynamicCommunityDetector& base) {
	this->baseAlgos.push_back(&base);
}

void DynamicEnsemble::setFinalAlgorithm(Clusterer& final) {
	this->finalAlgo = &final;
}

void DynamicEnsemble::setOverlapper(Overlapper& overlap) {
	this->overlapAlgo = &overlap;
}

} /* namespace NetworKit */
